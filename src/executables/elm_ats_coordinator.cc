
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
                         Amanzi::Comm_ptr_type comm ) :
    Coordinator(parameter_list, S, comm) {}

void ELM_ATSCoordinator::setup() {
  Teuchos::TimeMonitor monitor(*setup_timer_);
  Coordinator::setup();
}

void ELM_ATSCoordinator::initialize() {

  Coordinator::initialize();

  // get the intial timestep
  //if (!restart_) {
  //  double dt = get_dt(false);
  //  S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
  //}

  // visualization at IC
  // for testing
  visualize();
  checkpoint();

  // is this needed?
  Teuchos::TimeMonitor cycle_monitor(*cycle_timer_);

  // Make sure times are set up correctly
  AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT)
                         - S_->get_time(Amanzi::Tags::CURRENT)) < 1.e-4);
}

bool ELM_ATSCoordinator::advance(double dt) {

  // run model for single timestep
  // order is important
  // assign dt and advance NEXT
  S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
  S_->advance_time(Amanzi::Tags::NEXT, dt);
  // get time from tags
  double t_old = S_->get_time(Amanzi::Tags::CURRENT);
  double t_new = S_->get_time(Amanzi::Tags::NEXT);
  // check that dt and time tags align
  AMANZI_ASSERT(std::abs((t_new - t_old) - dt) < 1.e-4);

  // advance pks
  auto fail = pk_->AdvanceStep(t_old, t_new, false);
  if (!fail) fail |= !pk_->ValidStep();

  if (!fail) {
    // commit the state, copying NEXT --> CURRENT
    pk_->CommitStep(t_old, t_new, Amanzi::Tags::NEXT);

  } else {
    // Failed the timestep.
    // Potentially write out failed timestep for debugging
    for (const auto& vis : failed_visualization_) WriteVis(*vis, *S_);

    // copy from old time into new time to reset the timestep
    pk_->FailStep(t_old, t_new, Amanzi::Tags::NEXT);
    S_->set_time(Amanzi::Tags::NEXT, S_->get_time(Amanzi::Tags::CURRENT));
  }

  if (fail) {
    Errors::Message msg("ELM_ATSCoordinator: Coordinator advance failed.");
    Exceptions::amanzi_throw(msg);
  }

  // set CURRENT = NEXT
  // state advance_cycle - necessary??
  S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
  S_->advance_cycle();

  // make observations, vis, and checkpoints
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
  visualize();
  checkpoint(); // checkpoint with the new dt

  return false;

}

} // namespace ATS
