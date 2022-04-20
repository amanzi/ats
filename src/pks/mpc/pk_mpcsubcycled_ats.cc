/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Class for subcycling a slave step within a master step.
  Assumes that intermediate_time() can be used (i.e. this is not nestable?)

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#include "pk_mpcsubcycled_ats.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
PK_MPCSubcycled_ATS::PK_MPCSubcycled_ATS(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
  PK(pk_tree, global_list, S, soln),
  MPC<PK>(pk_tree, global_list, S, soln),
  subcycling(true) {

  init_(S);

  S_ = S;

  // Master PK is the PK whose time step size sets the size, the slave is subcycled.
  master_ = plist_->get<int>("master PK index", 0);
  slave_ = master_ == 1 ? 0 : 1;

  if (sub_pks_.size() != 2 || master_ > 1) {
    Errors::Message message("PK_MPCSubcycled_ATS: only MPCs with two sub-PKs can currently be subcycled");
    std::cout<<name_<<" sub pks size "<<sub_pks_.size() <<" master "<<master_<<"\n";
    Exceptions::amanzi_throw(message);
  }

  // min dt allowed in subcycling
  min_dt_ = plist_->get<double>("minimum subcycled relative dt", 1.e-5);
  subcycling = plist_->get<bool>("subcycling", true);
}
  

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double PK_MPCSubcycled_ATS::get_dt() {
  master_dt_ = sub_pks_[master_]->get_dt();
  slave_dt_ = sub_pks_[slave_]->get_dt();
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  if (subcycling) return master_dt_;
  else return slave_dt_;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool PK_MPCSubcycled_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;

  // advance the master PK using the full step size
  fail = sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  master_dt_ = t_new - t_old;
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  // --etc: unclear if state should be commited?
  sub_pks_[master_]->CommitStep(t_old, t_new, S_);

  // advance the slave, subcycling if needed
  S_->set_intermediate_time(t_old);
  bool done = false;

  double dt_next = slave_dt_;
  double dt_done = 0.;
  while (!done) {
    // do not overstep
    if (t_old + dt_done + dt_next > t_new) {
      dt_next = t_new - t_old - dt_done;
    }

    // set the intermediate time
    S_->set_intermediate_time(t_old + dt_done + dt_next);

    // take the step
    fail = sub_pks_[slave_]->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // if success, commit the state and increment to next intermediate
      // -- etc: unclear if state should be commited or not?
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_);
      dt_done += dt_next;
    }

    // check for done condition
    done = (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) || // finished the step
        (dt_next  < min_dt_); // failed
  }

  if (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) {
    // done, success
    // --etc: unclear if state should be commited or not?
    CommitStep(t_old, t_new, S_);
    return false;
  } else {
    return true;
  }
}

}  // namespace Amanzi

