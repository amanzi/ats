/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  MPC for subcycling one PK relative to another.
*/

#include "mpc_subcycled.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCSubcycled::MPCSubcycled(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
  PK(pk_tree, global_list, S, soln),
  MPC<PK>(pk_tree, global_list, S, soln),
  subcycling_(true)
{
  init_();

  // Master PK is the PK whose time step size sets the size, the subcycled is subcycled.
  subcycled_ = plist_->get<int>("subcycled PK index", 1);
  standard_ = subcycled_ == 1 ? 0 : 1;

  if (sub_pks_.size() != 2 || subcycled_ > 1) {
    Errors::Message message("MPCSubcycled: only MPCs with two sub-PKs can currently be subcycled");
    Exceptions::amanzi_throw(message);
  }

  // min dt allowed in subcycling
  min_dt_ = plist_->get<double>("minimum subcycled relative dt", 1.e-5);
  subcycling = plist_->get<bool>("subcycling", true);
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCSubcycled::get_dt()
{
  standard_dt_ = sub_pks_[standard_]->get_dt();
  subcycled_dt_ = sub_pks_[subcycled_]->get_dt();
  if (subcycled_dt_ > standard_dt_) subcycled_dt_ = standard_dt_;

  if (subcycling) return standard_dt_;
  else return subcycled_dt_;
}


// -----------------------------------------------------------------------------
// Set standard dt
// -----------------------------------------------------------------------------
void MPCSubcycled::set_dt(double dt) {
  standard_dt_ = dt;
  sub_pks_[standard_]->set_dt(dt);

  if (!subcycling_) {
    subcycled_dt_ = dt;
    sub_pks_[subcycled_]->set_dt(dt);
  }
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool MPCSubcycled::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  standard_dt_ = t_new - t_old;
  if (subcycled_dt_ > standard_dt_) subcycled_dt_ = standard_dt_;

  // advance the standard PK using the full step size
  if (standard_ == 0) {
    fail = sub_pks_[standard_]->AdvanceStep(t_old, t_new, reinit);
    if (fail) return fail;
  }

  // advance the subcycled, subcycling if needed
  S_->set_intermediate_time(t_old);
  bool done = false;

  double dt_next = subcycled_dt_;
  double dt_done = 0.;
  while (!done) {
    // do not overstep
    if (t_old + dt_done + dt_next > t_new) {
      dt_next = t_new - t_old - dt_done;
    }

    // set the intermediate time
    S_->set_intermediate_time(t_old + dt_done + dt_next);

    // take the step
    fail = sub_pks_[subcycled_]->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // if success, commit the state and increment to next intermediate
      // -- etc: unclear if state should be commited or not?
      sub_pks_[subcycled_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, tag_next_);
      dt_done += dt_next;
    }

    // check for done condition
    done = (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) || // finished the step
        (dt_next  < min_dt_); // failed
  }

  if (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) {
    return false;
  } else {
    return true;
  }
}

}  // namespace Amanzi

