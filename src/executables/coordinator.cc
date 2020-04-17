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
  coordinator_init_();
}

void
Coordinator::coordinator_init_() {
  coordinator_list_ = Teuchos::sublist(parameter_list_, "cycle driver");
  read_parameter_list_();

  // create the time step manager
  tsm_ = Teuchos::rcp(new Amanzi::TimeStepManager());

  // create the checkpoint, vis, and observation objects
  checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(Teuchos::sublist(parameter_list_, "checkpoint"), comm_));
  observations_ = Teuchos::rcp(new Amanzi::UnstructuredObservations(Teuchos::sublist(parameter_list_, "observations"), comm_));
  
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

void Coordinator::setup() {
  // Require data for timestep control
  S_->Require<double>("time", "", "time");
  S_->Require<double>("time", "next", "time");
  S_->Require<int>("cycle", "", "cycle");
  S_->Require<int>("cycle", "next", "cycle");
  S_->Require<double>("dt", "", "dt");

  // Setup the PKs, which sets all requirements
  pk_->Setup();

  // Setup the State, which ensures compatibility of evaluators and then
  // allocates memory.
  S_->Setup();
}

void Coordinator::initialize() {
  Teuchos::OSTab tab = vo_->getOSTab();

  // Initialize the state (initializes all dependent variables).
  S_->set_time("", t0_);
  S_->set_cycle("", cycle0_);
  S_->GetW<double>("dt", "", "dt") = 0.;

  pk_->Initialize();
  S_->Initialize();

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

void Coordinator::finalize() {
  // Force checkpoint at the end of simulation, and copy to checkpoint_final
  pk_->CalculateDiagnostics("");
  WriteCheckpoint(*checkpoint_, *S_, true);

  // flush observations to make sure they are saved
  observations_->Flush();
}


void Coordinator::read_parameter_list_() {
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


bool Coordinator::advance(double dt) {
  S_->set_time("next", S_->time("") + dt);
  S_->set_cycle("next", S_->cycle("") + 1);

  bool fail = pk_->AdvanceStep("", "next");
  fail |= !pk_->ValidStep("", "next");

  if (fail) {
    pk_->FailStep("", "next");
    
  } else {
    // commit the state
    pk_->CommitStep("", "next");
    S_->set_time("", S_->time("next"));
    S_->set_cycle("", S_->cycle("next"));
    
    // make observations, vis, and checkpoints
    observations_->MakeObservations(*S_);
    visualize();
    checkpoint(dt);

  }
  return fail;
}



// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void Coordinator::cycle_driver() {
  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  // start at time t = t0 and initialize the state.
  {
    Teuchos::TimeMonitor monitor(*setup_timer_);
    setup();   
    initialize();
  }

  // get the intial timestep -- note, this would have to be fixed for a true restart
  double dt = get_dt(false);

  // visualization at IC
  visualize();
  checkpoint(dt);

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

      fail = advance(dt);
      dt = get_dt(fail);
    } // while not finished
  }

  // finalizing simulation
  Teuchos::TimeMonitor::summarize(*vo_->os());
  finalize();

} // cycle driver


} // close namespace Amanzi
