/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "mpc_surface_subsurface_helpers.hh"

#include "mpc_permafrost_split_flux_columns_subcycled.hh"

#include "PK_Physical.hh"

namespace Amanzi {

MPCPermafrostSplitFluxColumnsSubcycled::MPCPermafrostSplitFluxColumnsSubcycled(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPCPermafrostSplitFluxColumns(FElist, plist, S, solution)
{
  subcycled_timestep_type_ = plist_->get<std::string>("subcycling timestep type","surface star timestep");
  subcycled_target_time_ = plist_->get<double>("subcycling timestep target",3600);
};

double MPCPermafrostSplitFluxColumnsSubcycled::get_dt()
{
  if (subcycled_timestep_type_ == "global minimum") {
    double dt_l = 1.e99;
    for (auto pk : sub_pks_) {
      dt_l = std::min(pk->get_dt(), dt_l);
    }
    double dt_g;
    S_next_->GetMesh(Keys::getDomain(p_primary_variable_star_))->get_comm()->MinAll(&dt_l, &dt_g, 1);
    return dt_g;

  } else if (subcycled_timestep_type_ == "global target") {
    double dt_l = 1.e99;
    for (auto pk : sub_pks_) {
      dt_l = std::min(pk->get_dt(), dt_l);
    }
    double dt_g;
    dt_l = std::max(dt_l,subcycled_target_time_);

    S_next_->GetMesh(Keys::getDomain(p_primary_variable_star_))->get_comm()->MinAll(&dt_l, &dt_g, 1);
    return dt_g;

  } else if (subcycled_timestep_type_ == "surface star timestep") {
    return sub_pks_[0]->get_dt();

  } else {
    Errors::Message msg;
    msg << "MPCPermafrostSplitFluxColumnsSubcycled: unknown \"subcycling timestep type\" : \""
        << subcycled_timestep_type_ << "\"";
    Exceptions::amanzi_throw(msg);
  }
}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCPermafrostSplitFluxColumnsSubcycled::AdvanceStep(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  int my_pid = S_next_->GetMesh("surface_star")->get_comm()->MyPID();
  // Advance the star system
  bool fail = false;
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Beginning timestepping on surface star system" << std::endl;

  // Advance the surface star system
  double t_inner = t_old;
  bool done_sf = false;

  S_inter_->set_time(t_old);

  int count_sf = 0;
  while (!done_sf) {
    double dt_inner = std::min(sub_pks_[0]->get_dt(), t_new - t_inner);

    *S_next_->GetScalarData("dt", "coordinator") = dt_inner;
    S_next_->set_time(t_inner + dt_inner);
    bool fail_inner = sub_pks_[0]->AdvanceStep(t_inner, t_inner+dt_inner, false);
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  step failed? " << fail_inner << std::endl;
    bool valid_inner = sub_pks_[0]->ValidStep();
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "  step valid? " << valid_inner << std::endl;
    }

    if (fail_inner || !valid_inner) {
      dt_inner = sub_pks_[0]->get_dt();
      S_next_->AssignDomain(*S_inter_, "surface_star");
      S_next_->set_time(S_->get_time(tag_inter_));
      S_next_->set_cycle(S_inter_->cycle());

      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

    } else {
      sub_pks_[0]->CommitStep(t_inner, t_inner + dt_inner, S_next_);
      t_inner += dt_inner;
      if (t_inner >= t_new - 1.e-10) {
        done_sf = true;
      }

      S_inter_->AssignDomain(*S_next_, "surface_star");
      S_inter_->set_time(S_->get_time(tag_next_));
      S_inter_->set_cycle(S_next_->cycle());
      dt_inner = sub_pks_[0]->get_dt();
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
      }

    if (dt_inner < 1.e-4) {
      Errors::Message msg;
      msg << "Surface_star crashing timestep in subcycling: dt = " << dt_inner;
      Exceptions::amanzi_throw(msg);
    }
  }
  S_inter_->set_time(t_old);

  // Copy star's new value into primary's old value
  CopyStarToPrimary(t_new - t_old);

  // Now advance the columns
  const auto& col_domain_set = *S_->GetDomainSet(domain_col_);
  auto ds_iter = col_domain_set.begin();
  for (int i=1; i!=sub_pks_.size(); ++i) {
    auto col_domain = *ds_iter;
    double t_inner = t_old;
    bool done = false;
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "Beginning timestepping on " << col_domain << std::endl;

    S_inter_->set_time(t_old);
    while (!done) {
      double dt_inner = std::min(sub_pks_[i]->get_dt(), t_new - t_inner);
      *S_next_->GetScalarData("dt", "coordinator") = dt_inner;
      S_next_->set_time(t_inner + dt_inner);
      bool fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner+dt_inner, false);
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  step failed? " << fail_inner << std::endl;
      bool valid_inner = sub_pks_[i]->ValidStep();
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << "  step valid? " << valid_inner << std::endl;

      // DEBUGGING
      std::cout << col_domain << " (" << my_pid << ") Step: " << t_inner/86400.0 << " (" << dt_inner/86400.
                << ") failed/!valid = " << fail_inner << "," << !valid_inner << std::endl;
      // END DEBUGGING
      }

      if (fail_inner || !valid_inner) {
        dt_inner = sub_pks_[i]->get_dt();
        S_next_->AssignDomain(*S_inter_, col_domain);
        S_next_->AssignDomain(*S_inter_, "surface_"+col_domain);
        S_next_->AssignDomain(*S_inter_, "snow_"+col_domain);
        S_next_->set_time(S_->get_time(tag_inter_));
        S_next_->set_cycle(S_inter_->cycle());

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

      } else {
        sub_pks_[i]->CommitStep(t_inner, t_inner + dt_inner, S_next_);
        t_inner += dt_inner;
        if (t_inner >= t_new - 1.e-10) done = true;

        S_inter_->AssignDomain(*S_next_, col_domain);
        S_inter_->AssignDomain(*S_next_, "surface_"+col_domain);
        S_inter_->AssignDomain(*S_next_, "snow_"+col_domain);
        S_inter_->set_time(S_->get_time(tag_next_));
        S_inter_->set_cycle(S_next_->cycle());
        dt_inner = sub_pks_[i]->get_dt();
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
      }

      if (dt_inner < 1.e-4) {
        Errors::Message msg;
        msg << "Column " << col_domain << " on PID " << my_pid << " crashing timestep in subcycling: dt = " << dt_inner;
        Exceptions::amanzi_throw(msg);
      }
    }
    ++ds_iter;
  }
  S_inter_->set_time(t_old);

  // Copy the primary into the star to advance
  CopyPrimaryToStar(S_next_.ptr(), S_next_.ptr());

  // this can never fail without error, because we always subcycle to ensure
  // the full timestep.
  return false;
}

bool MPCPermafrostSplitFluxColumnsSubcycled::ValidStep()
{
  // this is always valid, because the inner steps were valid
  return true;
}

void
MPCPermafrostSplitFluxColumnsSubcycled::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S)
{
  // the sub-PKs would call Commit a second time -- it was already called when
  // the inner step was successful and commited.  Calling it a second time is
  // bad because it messes up the nonlinear solver's inner solution history.
  return;
}


} // namespace
