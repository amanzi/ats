
#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "InputAnalysis.hh"

#include "Units.hh"
#include "CompositeVector.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "VisualizationDomainSet.hh"
#include "IO.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"


#include "pk_helpers.hh"
#include "elm_ats_coordinator.hh"

namespace ATS {


  ELM_ATSCoordinator::ELM_ATSCoordinator(Teuchos::ParameterList& parameter_list,
                                         Teuchos::RCP<Amanzi::State>& S,
                                         Amanzi::Comm_ptr_type comm )
  : Coordinator(parameter_list, S, comm) {}

void ELM_ATSCoordinator::setup() {
  Teuchos::TimeMonitor monitor(*setup_timer_);
  Coordinator::setup();
}

void ELM_ATSCoordinator::initialize() {

  Coordinator::initialize();

  // visualization at IC  - This will cause error in ::reinit()
  // for testing
  visualize();
  //checkpoint();

  // is this needed?
  Teuchos::TimeMonitor cycle_monitor(*cycle_timer_);

  // Make sure times are set up correctly
  AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT)
                         - S_->get_time(Amanzi::Tags::CURRENT)) < 1.e-4);
}

//
void ELM_ATSCoordinator::reinit(double start_time, bool visout) {
  //
  t0_ = start_time;

  S_->set_time(Amanzi::Tags::CURRENT, t0_);
  S_->set_time(Amanzi::Tags::NEXT, t0_);
  S_->set_cycle(cycle0_);

  pk_->CommitStep(t0_, t0_, Amanzi::Tags::NEXT);

  // commit the initial conditions
  pk_->CommitStep(S_->get_time(), S_->get_time(), Amanzi::Tags::NEXT);

  // visualization at reset ICs
  visualize(visout);

  // is this needed?
  Teuchos::TimeMonitor cycle_monitor(*cycle_timer_);

  // Make sure times are set up correctly
  AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT)
                         - S_->get_time(Amanzi::Tags::CURRENT)) < 1.e-4);
}

bool ELM_ATSCoordinator::advance(double dt, bool visout, bool chkout) {    
  // start and end times for timestep
  double t_end = S_->get_time() + dt;

  bool fail = false;
  while (S_->get_time() < t_end && dt > 0.0) {

    // run model for a duration of dt
    // order is important
    // advance NEXT time tag
    S_->advance_time(Amanzi::Tags::NEXT, dt);
    // get time from tags
    double t_old = S_->get_time(Amanzi::Tags::CURRENT);
    double t_new = S_->get_time(Amanzi::Tags::NEXT);
    // check that dt and time tags align
    AMANZI_ASSERT(std::abs((t_new - t_old) - dt) < 1.e-4);

    // advance pks
    fail = pk_->AdvanceStep(t_old, t_new, false);
    if (!fail) fail |= !pk_->ValidStep();

    if (fail) {

      // Failed the timestep.
      // Potentially write out failed timestep for debugging
      for (const auto& vis : failed_visualization_) WriteVis(*vis, *S_);

      // set as failed and revert timestamps
      pk_->FailStep(t_old, t_new, Amanzi::Tags::NEXT);
      // set NEXT = CURRENT to reset t_new
      S_->set_time(Amanzi::Tags::NEXT, S_->get_time(Amanzi::Tags::CURRENT));

    } else {

      // commit the state, copying NEXT --> CURRENT
      pk_->CommitStep(t_old, t_new, Amanzi::Tags::NEXT);

      // set CURRENT = NEXT
      S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
      // state advance_cycle - necessary??
      S_->advance_cycle();

      // make observations, vis, and checkpoints
      for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
      double t1 = get_end_time();
      if (t_end>=t1 && visout) {visualize(visout);}
      if (t_end>=t1 && chkout) {checkpoint(chkout);} // checkpoint with the new dt
    }

    // get new dt and assign to State
    dt = get_dt(fail);
    S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
  }

  return fail;
}

void ELM_ATSCoordinator::finalize()
{
  WriteStateStatistics(*S_, *vo_);
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());
  Coordinator::finalize();
}

} // namespace ATS
