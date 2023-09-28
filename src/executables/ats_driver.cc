/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include <Epetra_MpiComm.h>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "VerboseObject.hh"

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
#include "exceptions.hh"
#include "errors.hh"

#include "ats_driver.hh"

// won't run if DEBUG_MODE == false
// error: 'Exceptions' in namespace 'Amanzi' does not name a type
// 125   } catch (Amanzi::Exceptions::Amanzi_exception &e) {
#define DEBUG_MODE 1

namespace ATS {

// -----------------------------------------------------------------------------
// setup and initialize, then run until time >= duration_
// -----------------------------------------------------------------------------
void
ATSDriver::cycle_driver()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  //
  // setup phase
  // ----------------------------------------
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning setup stage..." << std::endl
               << std::flush;
  }
  setup();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("2: setup");
  }

  //
  // initialize phase
  // ----------------------------------------
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning initialize stage..." << std::endl
               << std::flush;
  }
  initialize();

  // get the intial timestep
  double dt = get_dt(false);
  S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt);

  // Write dependency graph, initial state
  S_->WriteDependencyGraph();
  WriteStateStatistics(*S_, *vo_);

  // checkpoint, observe, and vis at initial time
  visualize();
  observe();
  checkpoint();

  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("3: initialize");
  }

  //
  // timestepping phase
  // ----------------------------------------
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning timestepping stage..." << std::endl
               << std::flush;
  }

  // iterate process kernels
  //
  // Make sure times are set up correctly
  {
    Teuchos::TimeMonitor timer(*timers_.at("4: solve"));
    AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT) - S_->get_time(Amanzi::Tags::CURRENT)) <
                  1.e-4);
    AMANZI_ASSERT(std::abs(dt - S_->Get<double>("dt", Amanzi::Tags::DEFAULT)) < 1.e-4);

#if !DEBUG_MODE
    try {
#endif
      while (((t1_ < 0) || (S_->get_time() < t1_)) &&
             ((cycle1_ == -1) || (S_->get_cycle() <= cycle1_)) &&
             ((duration_ < 0) || (wallclock_timer_->totalElapsedTime(true) < duration)) &&
             (dt > 0.)) {
        if (vo_->os_OK(Teuchos::VERB_LOW)) {
          *vo_->os()
            << "================================================================================"
            << std::endl
            << std::endl
            << "Cycle = " << S_->get_cycle() << ",  Time [days] = " << std::setprecision(16)
            << S_->get_time() / (60 * 60 * 24) << ",  dt [days] = " << std::setprecision(16)
            << dt / (60 * 60 * 24) << std::endl
            << "--------------------------------------------------------------------------------"
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

          visualize();
          observe();
          checkpoint();
        }

        dt = get_dt(fail);
      } // while not finished

#if !DEBUG_MODE
    } catch (Errors::TimeStepCrash& e) {
      // write one more vis for help debugging
      S_->advance_cycle(Amanzi::Tags::NEXT);
      visualize(true); // force vis

      // flush observations to make sure they are saved
      for (const auto& obs : observations_) obs->Flush();

      // dump a post_mortem checkpoint file for debugging
      checkpoint_->set_filebasename("post_mortem");
      checkpoint_->Write(*S_, Amanzi::Checkpoint::WriteType::POST_MORTEM);
      throw e;
    }
#endif
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("4: solve");
  }


  // finalizing simulation
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning finalize stage..." << std::endl
               << std::flush;
  }
  finalize();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("5: finalize");
  }
} // cycle driver


// -----------------------------------------------------------------------------
// run simulation
// -----------------------------------------------------------------------------
int
ATSDriver::run()
{
  // run the simulation
  cycle_driver();
  return 0;
}


} // end namespace ATS
