/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A DomainSet coupler, couples a bunch of domains of the same structure.

------------------------------------------------------------------------- */

#include "mpc_domain_set.hh"

namespace Amanzi {

MPCDomainSet::MPCDomainSet(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution)
    : MPC<PK>(pk_tree, global_list, S, solution),
      PK(pk_tree, global_list, S, solution),
      subcycled_(false),
      subcycled_target_dt_(-1.)
{
  // grab the list of subpks
  auto subpks = this->plist_->template get<Teuchos::Array<std::string> >("PKs order");
  if (subpks.size() != 1) {
    Errors::Message msg;
    msg << "MPCDomainSet: \"" << name()
        << "\" expected exactly one entry in PKs order and it must be a domain set.";
    Exceptions::amanzi_throw(msg);
  }

  auto pk_name = subpks.back();
  subpks.pop_back();

  KeyTriple triple;
  bool is_ds = Keys::splitDomainSet(pk_name, triple);
  ds_name_ = std::get<0>(triple);
  if (!is_ds || !S->HasDomainSet(ds_name_)) {
    Errors::Message msg;
    msg << "MPCDomainSet: \"" << ds_name_ << "\" should be a domain-set of the form DOMAIN_*-PK_NAME";
    Exceptions::amanzi_throw(msg);
  }

  // add for the various sub-pks based on IDs
  auto ds = S->GetDomainSet(ds_name_);
  for (auto& subdomain : *ds) {
    subpks.push_back(Keys::getKey(subdomain, std::get<2>(triple)));
  }
  this->plist_->template set("PKs order", subpks);

  // construct the sub-PKs on COMM_SELF
  // FIXME: This should somehow get run on the DomainSet's entries Comms, not
  // on COMM_SELF! --etc
  MPC<PK>::init_(getCommSelf());

  // check whether we are subcycling
  subcycled_ = plist_->template get<bool>("subcycle subdomains", false);
  if (subcycled_) {
    subcycled_target_dt_ = plist_->template get<double>("subcycling target time step [s]");
    subcycled_min_dt_ = plist_->template get<double>("minimum subcycled time step [s]", 1.e-4);
  }
}


// must communicate dts since columns are serial
double MPCDomainSet::get_dt()
{
  double dt = 1.0e99;
  if (subcycled_) {
    dt = subcycled_target_dt_;
  } else {
    for (const auto& pk : sub_pks_) {
      dt = std::min<double>(dt,pk->get_dt());
    }

    double dt_local = dt;
    solution_->Comm()->MinAll(&dt_local, &dt, 1);
  }
  return dt;
}

// -----------------------------------------------------------------------------
// Set timestep for sub PKs
// -----------------------------------------------------------------------------
void MPCDomainSet::set_dt(double dt)
{
  if (subcycled_) {
    cycle_dt_ = dt;
  } else {
    for (const auto& pk : sub_pks_) {
      pk->set_dt(dt);
    }
  }
};


// -----------------------------------------------------------------------------
// Set tags for this and for subcycling
// -----------------------------------------------------------------------------
void MPCDomainSet::set_tags(const Tag& current, const Tag& next)
{
  if (subcycled_) {
    PK::set_tags(current, next);
    for (auto& pk : sub_pks_) {
      // create tags for subcycling
      Tag tag_subcycle_current(tag_current_.get()+" "+pk->name());
      Tag tag_subcycle_next(tag_next_.get()+" "+pk->name());
      pk->set_tags(tag_subcycle_current, tag_subcycle_next);
    }
  } else {
    MPC<PK>::set_tags(current, next);
  }
}



//-------------------------------------------------------------------------------------
// Advance the timestep
//-------------------------------------------------------------------------------------
bool
MPCDomainSet::AdvanceStep(double t_old, double t_new, bool reinit)
{
  if (subcycled_) return AdvanceStep_Subcycled_(t_old, t_new, reinit);
  else return AdvanceStep_Standard_(t_old, t_new, reinit);
}


//-------------------------------------------------------------------------------------
// Advance the timestep in the standard MPC way, but make sure to communicate failure
//-------------------------------------------------------------------------------------
bool
MPCDomainSet::AdvanceStep_Standard_(double t_old, double t_new, bool reinit)
{
  int nfailed = 0;
  for (const auto& pk : sub_pks_) {
    bool fail = pk->AdvanceStep(t_old, t_new, reinit);
    if (fail) {
      nfailed++;
      break;
    }
  }
  int nfailed_global(0);
  solution_->Comm()->SumAll(&nfailed, &nfailed_global, 1);
  if (nfailed_global) return true;
  return false;
}


//-------------------------------------------------------------------------------------
// Advance the timestep through subcyling
//-------------------------------------------------------------------------------------
bool
MPCDomainSet::AdvanceStep_Subcycled_(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  bool fail = false;
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Beginning subcycled timestepping." << std::endl;

  const auto& domain_set = *S_->GetDomainSet(ds_name_);
  int my_pid = solution_->Comm()->MyPID();

  int i = 0;
  for (const auto& domain : domain_set) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "Beginning subcyling on pk \"" << sub_pks_[i]->name() << "\"" << std::endl;

    bool done = false;
    double t_inner = t_old;

    Tag tag_subcycle_current(tag_current_.get()+" "+sub_pks_[i]->name());
    Tag tag_subcycle_next(tag_next_.get()+" "+sub_pks_[i]->name());

    S_->set_time(tag_subcycle_current, t_old);
    while (!done) {
      double dt_inner = std::min(sub_pks_[i]->get_dt(), t_new - t_inner);
      S_->Assign("dt", tag_subcycle_next, "coordinator", dt_inner);
      S_->set_time(tag_subcycle_next, t_inner + dt_inner);
      bool fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner+dt_inner, false);
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  step failed? " << fail_inner << std::endl;
      bool valid_inner = sub_pks_[i]->ValidStep();
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << "  step valid? " << valid_inner << std::endl;
      }

      if (fail_inner || !valid_inner) {
        sub_pks_[i]->FailStep(t_old, t_new, tag_subcycle_next);

        dt_inner = sub_pks_[i]->get_dt();
        S_->set_time(tag_subcycle_next, S_->get_time(tag_subcycle_current));

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

      } else {
        sub_pks_[i]->CommitStep(t_inner, t_inner + dt_inner, tag_subcycle_next);
        t_inner += dt_inner;
        if (std::abs(t_new - t_inner) < 1.e-10) done = true;

        S_->set_time(tag_subcycle_current, S_->get_time(tag_subcycle_next));
        S_->advance_cycle(tag_subcycle_current);

        dt_inner = sub_pks_[i]->get_dt();
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
      }

      if (dt_inner < subcycled_min_dt_) {
        Errors::Message msg;
        msg << "Subdomain " << domain << " on PID " << my_pid << " crashing timestep in subcycling: dt = " << dt_inner;
        Exceptions::amanzi_throw(msg);
      }
    }
    i++;
  }

  return false;
}


bool
MPCDomainSet::ValidStep()
{
  if (subcycled_) return true; // this was already checked in advance
  else {
    int valid_l = MPC<PK>::ValidStep() ? 1 : 0;
    int valid_g(-1);
    solution_->Comm()->MinAll(&valid_l, &valid_g, 1);
    return (valid_g == 1);
  }
}


void
MPCDomainSet::CommitStep(double t_old, double t_new, const Tag& tag)
{
  // In the case of subcycling, the sub-PKs would call Commit a second time --
  // it was already called when the inner step was successful and commited.
  // Calling it a second time is bad because it messes up the nonlinear
  // solver's inner solution history.
  if (subcycled_) return;
  else MPC<PK>::CommitStep(t_old, t_new, tag);
}

} // namespace Amanzi
