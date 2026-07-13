/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include <Epetra_MpiComm.h>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "VerboseObject.hh"

#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "TimeStepManager.hh"
#include "State.hh"
#include "PK.hh"
#include "IO.hh"

#include "exceptions.hh"
#include "errors.hh"

#include "time_advancer.hh"
#include "ats_driver.hh"

// won't run if DEBUG_MODE == false
// error: 'Exceptions' in namespace 'Amanzi' does not name a type
// 125   } catch (Amanzi::Exceptions::Amanzi_exception &e) {
#define DEBUG_MODE 1

namespace ATS {

void
ATSDriver::createTimeAdvancer_()
{
  // build a combined plist: done conditions + dt bounds from coordinator_list_,
  // vis/obs/checkpoint from plist_ top-level sublists
  auto ta_plist = Teuchos::rcp(new Teuchos::ParameterList(*coordinator_list_));
  if (plist_->isSublist("checkpoint"))
    ta_plist->sublist("checkpoint") = plist_->sublist("checkpoint");
  if (plist_->isSublist("observations"))
    ta_plist->sublist("observations") = plist_->sublist("observations");
  if (plist_->isSublist("visualization"))
    ta_plist->sublist("visualization") = plist_->sublist("visualization");
  if (plist_->isSublist("visualization failed"))
    ta_plist->sublist("visualization failed") = plist_->sublist("visualization failed");

  time_advancer_ = Teuchos::rcp(new TimeAdvancer(
    ta_plist, S_, pk_, tsm_,
    Amanzi::Tags::CURRENT, Amanzi::Tags::NEXT,
    vo_, wallclock_timer_));
}


void
ATSDriver::cycle_driver_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  //
  // setup phase
  // ----------------------------------------
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning setup stage..." << std::endl
               << std::flush;
  }
  parseParameterList();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("2a: parseParameterList");
  }

  // TimeAdvancer must be created after parseParameterList (tags/PKs set up)
  // but before setup() so observations can call Setup(S_) before State::Setup()
  createTimeAdvancer_();

  setup();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("2b: setup");
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

  S_->WriteDependencyGraph();
  WriteStateStatistics(*S_, *vo_);

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

  {
    Teuchos::TimeMonitor timer(*timers_.at("4: solve"));

#if !DEBUG_MODE
    try {
#endif
      time_advancer_->advance(t0_, t1_);
#if !DEBUG_MODE
    } catch (Errors::TimestepCrash& e) {
      // write one more vis for help debugging
      S_->advance_cycle(Amanzi::Tags::NEXT);
      time_advancer_->visualize(true);
      time_advancer_->observe();

      // dump a post_mortem checkpoint for debugging
      // TODO: expose post_mortem write through TimeAdvancer
      throw e;
    }
#endif
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("4: solve");
  }

  //
  // finalize phase
  // ----------------------------------------
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
}


int
ATSDriver::run()
{
  cycle_driver_();
  return 0;
}

} // end namespace ATS
