/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Base class for top-level ATS simulation drivers.
/*!

Driver is responsible for constructing meshes, State, and the PK hierarchy,
and for driving the PK lifecycle (parseParameterList, setup, initialize,
finalize).  Timestep advancement is delegated to TimeAdvancer.

Subclasses add a run() method (ATSDriver) or expose a per-step advance()
method (ELM_ATSDriver).

*/
#pragma once

#include <map>
#include <string>

#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_MpiComm.h"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "VerboseObject.hh"

namespace Amanzi {
namespace Utils {
class TimeStepManager;
} // namespace Utils
class State;
class TreeVector;
class PK;
} // namespace Amanzi

namespace ATS {

class TimeAdvancer;

class Driver {
 public:
  Driver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
         const Teuchos::RCP<Teuchos::Time>& wallclock_timer,
         const Teuchos::RCP<const Teuchos::Comm<int>>& teuchos_comm,
         const Amanzi::Comm_ptr_type& comm);

  virtual void parseParameterList();
  virtual void setup();
  virtual void initialize();
  virtual void finalize(bool checkpoint = true);
  void report_memory();

 protected:
  void initializeFromPlist_();
  void reportOneTimer_(const std::string& timer);

  Teuchos::RCP<Amanzi::PK> pk_;
  Teuchos::RCP<Amanzi::State> S_;
  Teuchos::RCP<Amanzi::TreeVector> soln_;
  Teuchos::RCP<Amanzi::Utils::TimeStepManager> tsm_;
  Teuchos::RCP<TimeAdvancer> time_advancer_;

  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;

  double t0_, t1_;
  int cycle0_;

  bool restart_;
  std::string restart_filename_;

  Amanzi::Comm_ptr_type comm_;
  Teuchos::RCP<const Teuchos::Comm<int>> teuchos_comm_;

  std::map<std::string, Teuchos::RCP<Teuchos::Time>> timers_;
  Teuchos::RCP<Teuchos::Time> wallclock_timer_;
  Teuchos::TimeMonitor wallclock_monitor_;

  Teuchos::RCP<Amanzi::VerboseObject> vo_;
};

} // namespace ATS
