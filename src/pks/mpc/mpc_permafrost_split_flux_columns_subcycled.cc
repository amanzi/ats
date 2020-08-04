/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"
#include "mpc_surface_subsurface_helpers.hh"

#include "mpc_permafrost_split_flux_columns_subcycled.hh"

#include "PK_Physical.hh"

namespace Amanzi {

MPCPermafrostSplitFluxColumnsSubcycled::MPCPermafrostSplitFluxColumnsSubcycled(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPCPermafrostSplitFluxColumns(FElist, plist, S, solution) {};


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
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  fail |= !sub_pks_[0]->ValidStep();
  if (fail) return fail;
  sub_pks_[0]->CommitStep(t_old, t_new, S_next_);

  // Copy star's new value into primary's old value
  CopyStarToPrimary(t_new - t_old);

  // Now advance the primary
  for (int i=1; i!=sub_pks_.size(); ++i) {
    auto col_domain = col_domains_[i-1];
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
        //S_next_->AssignDomain(*S_inter_, "surface_star");
        S_next_->set_time(S_inter_->time());
        S_next_->set_cycle(S_inter_->cycle());
        //*S_next_ = *S_inter_;

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;
        
      } else {
        sub_pks_[i]->CommitStep(t_inner, t_inner + dt_inner, S_next_);
        t_inner += dt_inner;
        if (t_inner >= t_new - 1.e-10) {
          done = true;
        }

        S_inter_->AssignDomain(*S_next_, col_domain);
        S_inter_->AssignDomain(*S_next_, "surface_"+col_domain);
        S_inter_->AssignDomain(*S_next_, "snow_"+col_domain);
        //        S_inter_->AssignDomain(*S_next_, "surface_star");
        S_inter_->set_time(S_next_->time());
        S_inter_->set_cycle(S_next_->cycle());
        // *S_inter_ = *S_next_;
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
  }
  S_inter_->set_time(t_old);

  // Copy the primary into the star to advance
  CopyPrimaryToStar(S_next_.ptr(), S_next_.ptr());

  return false;
}

bool MPCPermafrostSplitFluxColumnsSubcycled::ValidStep() 
{
  return true;
}
  


void MPCPermafrostSplitFluxColumnsSubcycled::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S)
{}

} // namespace
