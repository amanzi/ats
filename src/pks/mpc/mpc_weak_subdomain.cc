/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Weak MPC for subdomain model MPCs.
/*!

Weakly couples N PKs of the same type across a domain set.  Handles several
options in subcycling, subcommunicators, and other complexity associated with
this task.

Note that, unlike mpc_subcycled, this does not let a subset of PKs be
subcycled.  That seems appropriate as this means to couple PKs of the same type
-- why would you want to subcycle some columns but not all, or some watersheds
but not all?  That could be generalized, but it would be tricky to define on an
input spec.

*/


#include "mpc_weak_subdomain.hh"


namespace Amanzi {


MPCWeakSubdomain::MPCWeakSubdomain(Teuchos::ParameterList& FElist,
                                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& solution)
  : PK(FElist, plist, S, solution), MPC<PK>(FElist, plist, S, solution)
{
  init_();

  // check whether we are subcycling
  subcycled_ = plist_->template get<bool>("subcycle", false);
  if (subcycled_) {
    subcycled_target_dt_ = plist_->template get<double>("subcycling target time step [s]");
  }
};


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
MPCWeakSubdomain::get_dt()
{
  double dt = std::numeric_limits<double>::max();

  if (subcycled_) {
    dt = subcycled_target_dt_;
  } else {
    for (auto& pk : sub_pks_) dt = std::min(dt, pk->get_dt());
    double dt_local = dt;
    solution_->Comm()->MinAll(&dt_local, &dt, 1);
  }
  return dt;
}

// -----------------------------------------------------------------------------
// Set timestep for sub PKs
// -----------------------------------------------------------------------------
void
MPCWeakSubdomain::set_dt(double dt)
{
  if (subcycled_) {
    cycle_dt_ = dt;
  } else {
    for (auto& pk : sub_pks_) pk->set_dt(dt);
  }
};


// -----------------------------------------------------------------------------
// Set tags for this and for subcycling
// -----------------------------------------------------------------------------
void
MPCWeakSubdomain::set_tags(const Tag& current, const Tag& next)
{
  if (subcycled_) {
    PK::set_tags(current, next);

    const auto& ds = S_->GetDomainSet(ds_name_);
    int i = 0;
    for (const auto& subdomain : *ds) {
      // create tags for subcycling
      Tag tag_subcycle_current = get_ds_tag_current_(subdomain);
      Tag tag_subcycle_next = get_ds_tag_next_(subdomain);
      sub_pks_[i]->set_tags(tag_subcycle_current, tag_subcycle_next);
      ++i;
    }
  } else {
    MPC<PK>::set_tags(current, next);
  }
}


void
MPCWeakSubdomain::Setup()
{
  if (subcycled_) {
    const auto& ds = S_->GetDomainSet(ds_name_);
    for (const auto& subdomain : *ds) {
      Tag tag_subcycle_current = get_ds_tag_current_(subdomain);
      Tag tag_subcycle_next = get_ds_tag_next_(subdomain);

      S_->require_time(tag_subcycle_current);
      S_->require_time(tag_subcycle_next);
      S_->require_cycle(tag_subcycle_next);
      S_->Require<double>("dt", tag_subcycle_next, name());
    }
  }
  MPC<PK>::Setup();
}


void
MPCWeakSubdomain::Initialize()
{
  if (subcycled_) {
    const auto& ds = S_->GetDomainSet(ds_name_);
    for (const auto& subdomain : *ds) {
      Tag tag_subcycle_next = get_ds_tag_next_(subdomain);
      S_->GetRecordW("dt", tag_subcycle_next, name()).set_initialized();
    }
  }
  MPC<PK>::Initialize();
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool
MPCWeakSubdomain::AdvanceStep(double t_old, double t_new, bool reinit)
{
  if (subcycled_)
    return AdvanceStep_Subcycled_(t_old, t_new, reinit);
  else
    return AdvanceStep_Standard_(t_old, t_new, reinit);
}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool
MPCWeakSubdomain::AdvanceStep_Standard_(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  for (auto& pk : sub_pks_) {
    fail = pk->AdvanceStep(t_old, t_new, reinit);
    if (fail) break;
  }

  int sub_fail_i = fail ? 1 : 0;
  int fail_i;
  comm_->SumAll(&sub_fail_i, &fail_i, 1);
  return (bool)fail_i;
};


//-------------------------------------------------------------------------------------
// Advance the timestep through subcyling
//-------------------------------------------------------------------------------------
bool
MPCWeakSubdomain::AdvanceStep_Subcycled_(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  bool fail = false;
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Beginning subcycled timestepping." << std::endl;

  const auto& ds = *S_->GetDomainSet(ds_name_);
  int my_pid = solution_->Comm()->MyPID();

  int i = 0;
  int n_throw = 0;
  std::string throw_msg;
  for (const auto& subdomain : ds) {
    double dt_inner = -1;
    double t_inner = t_old;
    try { // must catch non-collective throws for TimeStepCrash
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "Beginning subcyling on pk \"" << sub_pks_[i]->name() << "\"" << std::endl;

      bool done = false;
      Tag tag_subcycle_current = get_ds_tag_current_(subdomain);
      Tag tag_subcycle_next = get_ds_tag_next_(subdomain);

      S_->set_time(tag_subcycle_current, t_old);
      while (!done) {
        dt_inner = std::min(sub_pks_[i]->get_dt(), t_new - t_inner);
        S_->Assign("dt", tag_subcycle_next, name(), dt_inner);
        S_->set_time(tag_subcycle_next, t_inner + dt_inner);
        bool fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner + dt_inner, false);
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
          S_->advance_cycle(tag_subcycle_next);

          dt_inner = sub_pks_[i]->get_dt();
          if (vo_->os_OK(Teuchos::VERB_EXTREME))
            *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
        }
      }
      i++;
    } catch (Errors::TimeStepCrash& e) {
      n_throw++;
      throw_msg = e.what();
      break;
    }
  }

  // check for any other ranks throwing and, if so, throw ourselves so that all procs throw
  int n_throw_g = 0;
  comm_->SumAll(&n_throw, &n_throw_g, 1);
  if (n_throw > 0) {
    // inject more information into the crash message
    Errors::TimeStepCrash msg;
    msg << throw_msg << "  on rank " << comm_->MyPID() << " of " << comm_->NumProc();
    Exceptions::amanzi_throw(msg);
  } else if (n_throw_g > 0) {
    Errors::TimeStepCrash msg;
    msg << "TimeStepCrash on another rank: nprocs failed = " << n_throw_g;
    Exceptions::amanzi_throw(msg);
  }
  return false;
}


void
MPCWeakSubdomain::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  if (S_->get_cycle() < 0 && tag_next == Tags::NEXT) {
    // initial commit, also do the substep commits
    if (subcycled_) {
      const auto& ds = S_->GetDomainSet(ds_name_);
      int i = 0;
      for (auto& domain : *ds) {
        auto l_tag_next = get_ds_tag_next_(domain);
        sub_pks_[i]->CommitStep(t_old, t_new, l_tag_next);
        ++i;
      }
    }
  }

  if (tag_next == tag_next_ && tag_next != Tags::NEXT) {
    // do not commit step in this case -- this is nested subcycling, which we
    // do not have a formal way of dealing with correctly.
    return;
  } else {
    for (const auto& pk : sub_pks_) { pk->CommitStep(t_old, t_new, tag_next); }
  }
}


void
MPCWeakSubdomain::init_()
{
  // grab the list of subpks
  auto subpks = plist_->get<Teuchos::Array<std::string>>("PKs order");
  if (subpks.size() != 1) {
    Errors::Message msg;
    msg << "MPCWeakSubdomain: \"PKs order\" should consist of a single domain set of sub-pks.";
    Exceptions::amanzi_throw(msg);
  }
  std::string subdomain_name = subpks[0];
  subpks.pop_back();

  KeyTriple subdomain_triple;
  bool is_ds = Keys::splitDomainSet(subdomain_name, subdomain_triple);
  if (!is_ds) {
    Errors::Message msg;
    msg << "MPCWeakSubdomain: subpk \"" << subdomain_name
        << "\" should be a domain-set PK of the form SUBDOMAIN_DOMAIN_NAME_*-NAME";
    Exceptions::amanzi_throw(msg);
  }

  // get the domain set and save the comm of the parent mesh for later
  ds_name_ = std::get<0>(subdomain_triple);
  const auto& ds = S_->GetDomainSet(ds_name_);
  comm_ = ds->getIndexingParent()->getComm();

  // -- create the lifted PKs
  PKFactory pk_factory;
  for (const auto& subdomain : *ds) {
    auto mesh = S_->GetMesh(subdomain);
    // create the solution vector, noting that these are created on the
    // communicator associated with the mesh of the subdomain, which may differ
    // from the coupler's comm.
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(mesh->getComm()));
    solution_->PushBack(pk_soln);

    // create the PK
    auto subpk = Keys::getKey(subdomain, std::get<2>(subdomain_triple));
    sub_pks_.emplace_back(pk_factory.CreatePK(subpk, pk_tree_, global_list_, S_, pk_soln));
  }
}

} // namespace Amanzi
