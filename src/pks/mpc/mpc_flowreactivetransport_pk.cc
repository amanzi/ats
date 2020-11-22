/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of Flow_PK with Transport_PK and Chemestry_PK
*/

#include "mpc_flowreactivetransport_pk.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
FlowReactiveTransport_PK_ATS::FlowReactiveTransport_PK_ATS(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln) :
    PK(pk_tree, global_list, S, soln),
    PK_MPCSubcycled_ATS(pk_tree, global_list, S, soln) { 

  name_ = pk_tree.name();
  const char* result = name_.data();

  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(name_,"->"); 
  if (res.end() - name_.end() != 0) boost::algorithm::erase_head(name_,  res.end() - name_.begin());

  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(global_list, "PKs", true);
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(pks_list, name_, true);

  vo_ = Teuchos::null;
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = pk_list -> sublist("verbose object");
  vo_ =  Teuchos::rcp(new VerboseObject("Flow&TransportPK", vlist)); 

  flow_timer_ = Teuchos::TimeMonitor::getNewCounter("flow");
  rt_timer_ = Teuchos::TimeMonitor::getNewCounter("reactive tranport");
  
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double FlowReactiveTransport_PK_ATS::get_dt() {
  double dt = Amanzi::PK_MPCSubcycled_ATS::get_dt();
  return dt;
}


// -----------------------------------------------------------------------------
// Set master dt
// -----------------------------------------------------------------------------
void FlowReactiveTransport_PK_ATS::set_dt(double dt) {
  master_dt_ = dt;
  sub_pks_[master_]->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Make necessary operatios by the end of the time steps.
// -----------------------------------------------------------------------------
void FlowReactiveTransport_PK_ATS::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  sub_pks_[slave_]->CommitStep(t_old, t_new, S);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool FlowReactiveTransport_PK_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;

  // advance the master PK using the full step size

  Teuchos::OSTab tab = vo_->getOSTab();

  Teuchos::RCP<Teuchos::TimeMonitor> local_flow_monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*flow_timer_));
  
  fail = sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  fail |= !sub_pks_[master_]->ValidStep();

  local_flow_monitor = Teuchos::null;
  
  if (fail) {

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      *vo_->os()<<"Master step is failed\n";
    }
    return fail;
  }

  //return fail;
  master_dt_ = t_new - t_old;
  sub_pks_[master_]->CommitStep(t_old, t_new, S_next_);


  Teuchos::RCP<const Field> field_tmp = S_->GetFieldCopy("mass_flux", "next_timestep");
  Key copy_owner = field_tmp->owner();
  Teuchos::RCP<Epetra_MultiVector> flux_copy = S_->GetFieldCopyData("mass_flux", "next_timestep", copy_owner)->ViewComponent("face", true);
  *flux_copy = *S_next_->GetFieldData("mass_flux")->ViewComponent("face", true);
 
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os()<<"Master step is successful\n";


  slave_dt_ = sub_pks_[slave_]->get_dt(); 
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;
  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) *vo_->os()<<"Slave dt="<<slave_dt_<<"\n";


  // advance the slave, subcycling if needed
  S_->set_initial_time(t_old);
  S_->set_intermediate_time(t_old);
  S_next_->set_intermediate_time(t_old);
  S_next_->set_final_time(t_new);
  bool done = false;

  double dt_next = slave_dt_;
  double dt_done = 0.;
  int ncycles = 0;

  local_flow_monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*rt_timer_));
  
  while (!done) {
    // do not overstep
    if (t_old + dt_done + dt_next > t_new) {
      dt_next = t_new - t_old - dt_done;
    }

    // take the step
    fail = sub_pks_[slave_]->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);
    ncycles ++;
    
    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // if success, commit the state and increment to next intermediate
      // -- etc: unclear if state should be commited or not?
      // set the intermediate time
      S_ -> set_intermediate_time(t_old + dt_done + dt_next);
      //S_next_ -> set_intermediate_time(t_old + dt_done + dt_next);
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_);
      //sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
      dt_done += dt_next;
    }

    // check for done condition
    done = (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) || // finished the step
        (dt_next  < min_dt_); // failed
  }


  if (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) {
    // done, success
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os()<<"Slave step is successful (use "<<ncycles<<" subcycles.)"<<"\n";
    return false;
  } else {
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os()<<"Slave step is failed\n";
    return true;
  }

  local_flow_monitor = Teuchos::null;
}

}  // namespace Amanzi

