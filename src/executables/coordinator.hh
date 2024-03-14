/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Simulation controller intended to be used as base class for top-level driver
#ifndef ATS_COORDINATOR_HH_
#define ATS_COORDINATOR_HH_

#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "VerboseObject.hh"

namespace Amanzi {
class TimeStepManager;
class Visualization;
class Checkpoint;
class State;
class TreeVector;
class PK;
class PK_ATS;
class UnstructuredObservations;
}; // namespace Amanzi


namespace ATS {

class Coordinator {
 public:
  Coordinator(const Teuchos::RCP<Teuchos::ParameterList>& plist,
              const Teuchos::RCP<Teuchos::Time>& wallclock_timer,
              const Teuchos::RCP<const Teuchos::Comm<int>>& teuchos_comm,
              const Amanzi::Comm_ptr_type& comm);

  // PK methods
  void setup();
  void initialize();
  void finalize();
  void report_memory();

  bool advance();
  bool visualize(bool force = false);
  void observe();
  bool checkpoint(bool force = false);
  double get_dt(bool after_fail = false);

 protected:
  void initializeFromPlist_();
  void reportOneTimer_(const std::string& timer);

  // PK container and factory
  Teuchos::RCP<Amanzi::PK> pk_;

  // states
  Teuchos::RCP<Amanzi::State> S_;
  Teuchos::RCP<Amanzi::TreeVector> soln_;

  // time step manager
  Teuchos::RCP<Amanzi::TimeStepManager> tsm_;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;

  double t0_, t1_;
  double max_dt_, min_dt_;
  int cycle0_, cycle1_;

  // Epetra communicator
  Amanzi::Comm_ptr_type comm_;

  // vis and checkpointing
  std::vector<Teuchos::RCP<Amanzi::Visualization>> visualization_;
  std::vector<Teuchos::RCP<Amanzi::Visualization>> failed_visualization_;
  Teuchos::RCP<Amanzi::Checkpoint> checkpoint_;
  bool restart_;
  std::string restart_filename_;

  // observations
  std::vector<Teuchos::RCP<Amanzi::UnstructuredObservations>> observations_;

  // Teuchos Communicator for timers... will go away in tpetra
  Teuchos::RCP<const Teuchos::Comm<int>> teuchos_comm_;

  // timers
  std::map<std::string, Teuchos::RCP<Teuchos::Time>> timers_;
  Teuchos::RCP<Teuchos::Time> wallclock_timer_;
  Teuchos::TimeMonitor wallclock_monitor_;
  double duration_;
  bool subcycled_ts_;

  // fancy OS
  Teuchos::RCP<Amanzi::VerboseObject> vo_;
};


} // namespace ATS

#endif
