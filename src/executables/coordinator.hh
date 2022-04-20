/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Simulation controller and top-level driver

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

In the `"cycle driver`" sublist, the user specifies global control of the
simulation, including starting and ending times and restart options.

.. _coordinator-spec:
.. admonition:: coordinator-spec

    * `"start time`" ``[double]`` **0.** Specifies the start of time in model time.
    * `"start time units`" ``[string]`` **"s"** One of "s", "d", or "yr"

    ONE OF

    * `"end time`" ``[double]`` Specifies the end of the simulation in model time.
    * `"end time units`" ``[string]`` **"s"** One of `"s`", `"d`", or `"yr`"

    OR

    * `"end cycle`" ``[int]`` **optional** If provided, specifies the end of the
      simulation in timestep cycles.

      END

    * `"subcycled timestep`" ``[bool]`` **false**  If true, this coordinator creates
      a third State object to store intermediate solutions, allowing for failed
      steps.
    * `"restart from checkpoint file`" ``[string]`` **optional** If provided,
      specifies a path to the checkpoint file to continue a stopped simulation.
    * `"wallclock duration [hrs]`" ``[double]`` **optional** After this time, the
      simulation will checkpoint and end.
    * `"required times`" ``[io-event-spec]`` **optional** An IOEvent_ spec that
      sets a collection of times/cycles at which the simulation is guaranteed to
      hit exactly.  This is useful for situations such as where data is provided at
      a regular interval, and interpolation error related to that data is to be
      minimized.
    * `"PK tree`" ``[pk-typed-spec-list]`` List of length one, the top level
      PK_ spec.

Note: Either `"end cycle`" or `"end time`" are required, and if
both are present, the simulation will stop with whichever arrives
first.  An `"end cycle`" is commonly used to ensure that, in the case
of a time step crash, we do not continue on forever spewing output.

Example:

.. code-block:: xml

   <ParameterList name="cycle driver">
     <Parameter  name="end cycle" type="int" value="6000"/>
     <Parameter  name="start time" type="double" value="0."/>
     <Parameter  name="start time units" type="string" value="s"/>
     <Parameter  name="end time" type="double" value="1"/>
     <Parameter  name="end time units" type="string" value="yr"/>
     <ParameterList name="required times">
       <Parameter name="start period stop" type="Array(double)" value="{0,-1,86400}" />
     </ParameterList>
     <ParameterList name="PK tree">
       <ParameterList name="my richards pk">
         <Parameter name="PK type" type="string" value="richards" />
       </ParameterList>
     </ParameterList>
   </ParameterList>

*/

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
};


namespace ATS {

class Coordinator {

public:
  Coordinator(Teuchos::ParameterList& parameter_list,
              Teuchos::RCP<Amanzi::State>& S,
              Amanzi::Comm_ptr_type comm);
              //              Amanzi::ObservationData& output_observations);

  // PK methods
  void setup();
  double initialize(); // returns dt if a restarted run
  void finalize();
  void report_memory();
  bool advance(double t_old, double t_new, double& dt_next);
  void visualize(bool force=false);
  void checkpoint(double dt, bool force=false);
  double get_dt(bool after_fail=false);
  Teuchos::RCP<Amanzi::State> get_next_state() { return S_next_; }

  // one stop shopping
  void cycle_driver();

private:
  void coordinator_init();
  void read_parameter_list();

  // PK container and factory
  Teuchos::RCP<Amanzi::PK> pk_;

  // states
  Teuchos::RCP<Amanzi::State> S_;
  Teuchos::RCP<Amanzi::State> S_inter_;
  Teuchos::RCP<Amanzi::State> S_next_;
  Teuchos::RCP<Amanzi::TreeVector> soln_;

  // time step manager
  Teuchos::RCP<Amanzi::TimeStepManager> tsm_;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> parameter_list_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;

  double t0_, t1_;
  double max_dt_, min_dt_;
  int cycle0_, cycle1_;

  // Epetra communicator
  Amanzi::Comm_ptr_type comm_;

  // vis and checkpointing
  std::vector<Teuchos::RCP<Amanzi::Visualization> > visualization_;
  std::vector<Teuchos::RCP<Amanzi::Visualization> > failed_visualization_;
  Teuchos::RCP<Amanzi::Checkpoint> checkpoint_;
  bool restart_;
  std::string restart_filename_;

  // observations
  std::vector<Teuchos::RCP<Amanzi::UnstructuredObservations>> observations_;

  // timers
  Teuchos::RCP<Teuchos::Time> setup_timer_;
  Teuchos::RCP<Teuchos::Time> cycle_timer_;
  Teuchos::RCP<Teuchos::Time> timer_;
  double duration_;
  bool subcycled_ts_;

  // fancy OS
  Teuchos::RCP<Amanzi::VerboseObject> vo_;
};



} // close namespace ATS

#endif
