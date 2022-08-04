//duration_
//setup_timer_
//restart_
//S_
//cycle_timer_
//t1_
//cycle1_
//timer_

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

#include "ats_driver.hh"



namespace ATS {

ATSDriver::ATSDriver(Teuchos::ParameterList& parameter_list,
                 Teuchos::RCP<Amanzi::State>& S,
                 Amanzi::Comm_ptr_type comm)
    : Coordinator(parameter_list, S, comm) {}



// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void ATSDriver::cycle_driver() {
  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  // start at time t = t0 and initialize the state.
  {
    Teuchos::TimeMonitor monitor(*setup_timer_);
    setup();
    initialize();
  }

  // get the intial timestep
  double dt = get_dt(false);
  if (!restart_) {
    S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
  }

  // visualization at IC
  visualize();
  checkpoint();

  // iterate process kernels
  //
  // Make sure times are set up correctly
  AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT)
                         - S_->get_time(Amanzi::Tags::CURRENT)) < 1.e-4);
  {
    Teuchos::TimeMonitor cycle_monitor(*cycle_timer_);
    double dt = S_->Get<double>("dt", Amanzi::Tags::DEFAULT);
#if !DEBUG_MODE
  try {
#endif

    while (((t1_ < 0) || (S_->get_time() < t1_)) &&
           ((cycle1_ == -1) || (S_->get_cycle() <= cycle1_)) &&
           ((duration_ < 0) || (timer_->totalElapsedTime(true) < duration)) &&
           (dt > 0.)) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "======================================================================"
                  << std::endl << std::endl;
        *vo_->os() << "Cycle = " << S_->get_cycle();
        *vo_->os() << ",  Time [days] = "<< std::setprecision(16) << S_->get_time() / (60*60*24);
        *vo_->os() << ",  dt [days] = " << std::setprecision(16) << dt / (60*60*24)  << std::endl;
        *vo_->os() << "----------------------------------------------------------------------"
                  << std::endl;
      }

      S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);
      S_->advance_time(Amanzi::Tags::NEXT, dt);
      bool fail = advance();

      if (fail) {
        // reset t_new
        S_->set_time(Amanzi::Tags::NEXT, S_->get_time(Amanzi::Tags::CURRENT));
      } else {
        S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
        S_->advance_cycle();

        // make observations, vis, and checkpoints
        for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
        visualize();
        checkpoint(); // checkpoint with the new dt
      }

      dt = get_dt(fail);
    } // while not finished

#if !DEBUG_MODE
  } catch (Amanzi::Exceptions::Amanzi_exception &e) {
    // write one more vis for help debugging
    S_->advance_cycle(Amanzi::Tags::NEXT);
    visualize(true); // force vis

    // flush observations to make sure they are saved
    for (const auto& obs : observations_) obs->Flush();

    // catch errors to dump two checkpoints -- one as a "last good" checkpoint
    // and one as a "debugging data" checkpoint.
    checkpoint_->set_filebasename("last_good_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), comm_, *S_);
    checkpoint_->set_filebasename("error_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), comm_, *S_);
    throw e;
  }
#endif
  }

  // finalizing simulation
  WriteStateStatistics(*S_, *vo_);
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());

  finalize();
} // cycle driver


} // end namespace ATS
