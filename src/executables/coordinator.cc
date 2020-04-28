/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the Coordinator.  Coordinator is basically just a class to hold
the cycle driver, which runs the overall, top level timestep loop.  It
instantiates states, ensures they are initialized, and runs the timestep loop
including Vis and restart/checkpoint dumps.  It contains one and only one PK
-- most likely this PK is an MPC of some type -- to do the actual work.
------------------------------------------------------------------------- */

#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "Units.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"
#include "EvaluatorPrimary.hh"

#include "coordinator.hh"

#define DEBUG_MODE 1

namespace ATS {

Coordinator::Coordinator(const Teuchos::RCP<Teuchos::ParameterList>& parameter_list,
                         const Teuchos::RCP<Amanzi::State>& S,
                         const Amanzi::Comm_ptr_type& comm ) :
    parameter_list_(parameter_list),
    S_(S),
    comm_(comm),
    restart_(false) {

  // create and start the global timer
  timer_ = Teuchos::rcp(new Teuchos::Time("wallclock_monitor",true));
  setup_timer_ = Teuchos::TimeMonitor::getNewCounter("setup");
  cycle_timer_ = Teuchos::TimeMonitor::getNewCounter("cycle");

  vo_ = Teuchos::rcp(new Amanzi::VerboseObject("Coordinator", *parameter_list_));
  Init_();
}

void
Coordinator::Init_() {
  coordinator_list_ = Teuchos::sublist(parameter_list_, "cycle driver");
  ReadParameterList_();

  // create the time step manager
  tsm_ = Teuchos::rcp(new Amanzi::TimeStepManager());

  // create the checkpoint, vis, and observation objects
  checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(Teuchos::sublist(parameter_list_, "checkpoint"), comm_));
  observations_ = Teuchos::rcp(new Amanzi::UnstructuredObservations(Teuchos::sublist(parameter_list_, "observations"), comm_));
  auto vis_list = Teuchos::sublist(parameter_list_,"visualization");
  for (auto& entry : *vis_list) {
    std::string domain_name = entry.first;

    if (S_->HasMesh(domain_name)) {
      // visualize standard domain
      auto mesh_p = S_->GetMesh(domain_name);
      auto sublist_p = Teuchos::sublist(vis_list, domain_name);

      if (S_->HasMesh(domain_name+"_3d") && sublist_p->get<bool>("visualize on 3D mesh", true))
        mesh_p = S_->GetMesh(domain_name+"_3d");
      
      // vis successful timesteps
      auto vis = Teuchos::rcp(new Amanzi::Visualization(sublist_p, mesh_p));
      visualization_.push_back(vis);
    }
  }
  
  // create the top level PK
  auto pks_list = Teuchos::sublist(parameter_list_, "PKs");
  auto pk_tree_list = Teuchos::sublist(coordinator_list_, "PK tree");
  if (pk_tree_list->numParams() != 1) {
    Errors::Message message("CycleDriver: PK tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  auto pk_item = pk_tree_list->begin();
  auto pk_name = pk_tree_list->name(pk_item);
  
  // create the pk
  Amanzi::PKFactory pk_factory;
  pk_ = pk_factory.CreatePK(pk_name, pk_tree_list, parameter_list_, S_);

}

void Coordinator::Setup() {
  // Require data for timestep control
  S_->Require<double>("time", "", "time");
  S_->Require<double>("time", "next", "time");
  S_->Require<int>("cycle", "", "cycle");
  S_->Require<int>("cycle", "next", "cycle");
  S_->Require<double>("dt", "", "dt");

  // get primary variable evaluators for time
  Teuchos::ParameterList time_eval_plist("time");
  time_eval_plist.set("evaluator type", "primary variable double");
  S_->FEList().set("time", time_eval_plist);

  S_->RequireEvaluator("time", "");
  auto t0_eval = S_->GetEvaluatorPtr("time", "");
  t0_eval_ = Teuchos::rcp_dynamic_cast<Amanzi::EvaluatorPrimary<double>>(t0_eval);
  AMANZI_ASSERT(t0_eval_.get());

  S_->RequireEvaluator("time", "next");
  auto t1_eval = S_->GetEvaluatorPtr("time", "next");
  t1_eval_ = Teuchos::rcp_dynamic_cast<Amanzi::EvaluatorPrimary<double>>(t1_eval);
  AMANZI_ASSERT(t1_eval_.get());

  // Setup the PKs, which sets all requirements
  pk_->Setup();

  // Setup the State, which ensures compatibility of evaluators and then
  // allocates memory.
  S_->Setup();
}

void Coordinator::Initialize() {
  Teuchos::OSTab tab = vo_->getOSTab();

  // Initialize the state (initializes all dependent variables).
  S_->set_time("", t0_);
  t0_eval_->SetChanged();
  
  S_->set_cycle("", cycle0_);
  S_->GetW<double>("dt", "", "dt") = 0.;

  S_->Initialize();
  pk_->Initialize();

  // set up the TSM
  // -- register visualization times
  for (const auto& vis : visualization_) {
    vis->RegisterWithTimeStepManager(*tsm_);
  }

  // -- register checkpoint times
  checkpoint_->RegisterWithTimeStepManager(*tsm_);

  // -- register observation times
  observations_->RegisterWithTimeStepManager(*tsm_);

  // -- register the final time
  tsm_->RegisterTimeEvent(t1_);

  // -- register any intermediate requested times
  if (coordinator_list_->isSublist("required times")) {
    Amanzi::IOEvent pause_times(Teuchos::sublist(coordinator_list_, "required times"));
    pause_times.RegisterWithTimeStepManager(*tsm_);
  }
}

void Coordinator::Finalize() {
  // Force checkpoint at the end of simulation, and copy to checkpoint_final
  pk_->CalculateDiagnostics("");
  WriteCheckpoint(*checkpoint_, *S_, true);

  // flush observations to make sure they are saved
  observations_->Flush();
}


void Coordinator::ReadParameterList_() {
  Amanzi::Utils::Units units;
  t0_ = Amanzi::Utils::readValueAndUnits(units, *coordinator_list_, "start time", "s");
  t1_ = Amanzi::Utils::readValueAndUnits(units, *coordinator_list_, "end time", "s");
  max_dt_ = Amanzi::Utils::readValueAndUnits(units, *coordinator_list_,
          "max time step size", "s", 1.0e99);
  min_dt_ = Amanzi::Utils::readValueAndUnits(units, *coordinator_list_,
          "min time step size", "s", 1.0e-99);

  duration_ = Amanzi::Utils::readValueAndUnits(units, *coordinator_list_,
          "wallclock duration", "h", 24.0);

  cycle0_ = coordinator_list_->get<int>("start cycle",0);
  cycle1_ = coordinator_list_->get<int>("end cycle",-1);

  // restart control
  restart_ = coordinator_list_->isParameter("restart from checkpoint file");
  if (restart_) {
    restart_filename_ = coordinator_list_->get<std::string>("restart from checkpoint file");
  }
}


// -----------------------------------------------------------------------------
// acquire the chosen timestep size
// -----------------------------------------------------------------------------
double Coordinator::get_dt(bool after_fail) {
  // get the physical step size
  double dt = pk_->get_dt();

  // check if the step size has gotten too small
  if (dt < min_dt_) {
    Errors::Message message("Coordinator: error, timestep too small");
    Exceptions::amanzi_throw(message);
  }

  // cap the max step size
  if (dt > max_dt_) {
    dt = max_dt_;
  }

  // ask the step manager if this step size is ok
  dt = tsm_->TimeStep(S_->time("next"), dt, after_fail);
  return dt;
}


bool Coordinator::Advance(double dt) {
  S_->set_time("next", S_->time("") + dt);
  t1_eval_->SetChanged();  
  S_->set_cycle("next", S_->cycle("") + 1);

  bool fail = pk_->AdvanceStep("", "next");
  fail |= !pk_->ValidStep("", "next");

  if (fail) {
    pk_->FailStep("", "next");
    
  } else {
    // commit the state
    pk_->CommitStep("", "next");
    double time = S_->time("next");
    int cycle = S_->cycle("next");
    S_->set_time("", time);
    t0_eval_->SetChanged();  
    S_->set_cycle("", cycle);
    
    // make observations, vis, and checkpoints
    observations_->MakeObservations(*S_);
    for (const auto& v : visualization_) {
      if (v->DumpRequested(cycle, time)) WriteVis(*v, *S_);
    }
    if (checkpoint_->DumpRequested(cycle, time)) WriteCheckpoint(*checkpoint_, *S_);

  }
  return fail;
}



// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void Coordinator::CycleDriver() {
  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  // start at time t = t0 and initialize the state.
  {
    Teuchos::TimeMonitor monitor(*setup_timer_);
    Setup();   
    Initialize();
  }

  // get the intial timestep -- note, this would have to be fixed for a true restart
  double dt = get_dt(false);
  double time = S_->time("");
  int cycle = S_->cycle("");

  // visualization at IC
  for (const auto& v : visualization_) {
    if (v->DumpRequested(cycle, time)) WriteVis(*v, *S_);
  }
  if (checkpoint_->DumpRequested(cycle, time)) WriteCheckpoint(*checkpoint_, *S_);
  
  // iterate process kernels
  {
    Teuchos::TimeMonitor cycle_monitor(*cycle_timer_);

    bool fail = false;
    while ((S_->time() < t1_) &&
           ((cycle1_ == -1) || (S_->cycle() <= cycle1_)) &&
           (duration_ < 0 || timer_->totalElapsedTime(true) < duration) &&
           dt > 0.) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "======================================================================"
                  << std::endl << std::endl;
        *vo_->os() << "Cycle = " << S_->cycle();
        *vo_->os() << std::setprecision(15) << ",  Time [days] = "<< S_->time() / (60*60*24);
        *vo_->os() << ",  dt [days] = " << dt / (60*60*24)  << std::endl;
        *vo_->os() << "----------------------------------------------------------------------"
                  << std::endl;
      }

      S_->GetW<double>("dt", "", "dt") = dt;

      fail = Advance(dt);
      dt = get_dt(fail);
    } // while not finished
  }

  // finalizing simulation
  Teuchos::TimeMonitor::summarize(*vo_->os());
  Finalize();

} // cycle driver


} // close namespace Amanzi
