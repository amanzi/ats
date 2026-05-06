/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Drives a PK through a sequence of timesteps from t_start to t_end.
/*!

TimeAdvancer owns the while(!done) timestep loop, including vis, observations,
checkpointing, and deformable mesh recovery on failed steps.  It is used by
standalone drivers (ATSDriver), externally-coupled drivers (ELM_ATSDriver), and
subcycled MPCs (MPCSubcycled).

The single public entry point is:

   bool advance(double t_start, double t_end);

which drives the PK from t_start to t_end, subcycling internally as needed.

.. _time-advancer-spec:
.. admonition:: time-advancer-spec

   * `"end cycle`" ``[int]`` **-1** If >= 0, stop when cycle reaches this value.
   * `"wallclock duration [hrs]`" ``[double]`` **-1** If > 0, stop after this many
     wallclock hours.
   * `"max timestep size [s]`" ``[double]`` **1e99** Cap on dt.
   * `"min timestep size [s]`" ``[double]`` **1e-12** If dt falls below this, throw.
   * `"subcycled timestep`" ``[bool]`` **false** If true, limit dt to the PK's
     suggested dt even after TSM expansion.
   * `"checkpoint`" ``[checkpoint-spec]`` **optional** Checkpoint output spec.
   * `"observations`" ``[observation-spec-list]`` **optional** Observation specs.
   * `"visualization`" ``[visualization-spec-list]`` **optional** Vis specs.
   * `"visualization failed`" ``[visualization-spec-list]`` **optional** Vis on
     failed steps.

*/

#pragma once

#include <map>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"

#include "Key.hh"
#include "Tag.hh"
#include "VerboseObject.hh"

namespace Amanzi {
class PK;
class State;
class Visualization;
class Checkpoint;
class UnstructuredObservations;
namespace Utils {
class TimeStepManager;
} // namespace Utils
} // namespace Amanzi

namespace ATS {

class TimeAdvancer {
 public:
  TimeAdvancer(const Teuchos::RCP<Teuchos::ParameterList>& plist,
               const Teuchos::RCP<Amanzi::State>& S,
               const Teuchos::RCP<Amanzi::PK>& pk,
               const Teuchos::RCP<Amanzi::Utils::TimeStepManager>& tsm,
               const Amanzi::Tag& tag_current,
               const Amanzi::Tag& tag_next,
               const Teuchos::RCP<Amanzi::VerboseObject>& vo,
               const Teuchos::RCP<Teuchos::Time>& wallclock_timer = Teuchos::null);

  // Must be called during Driver::setup(), before State::Setup(), to allow
  // observations to require fields on State.
  void setup();

  // Must be called during Driver::initialize(), after State is fully
  // initialized, to dump IC output.
  void initialize();

  // Must be called during Driver::finalize() to flush observations and
  // optionally write a final checkpoint.
  void finalize(bool checkpoint = true);

  // Advance the PK from t_start to t_end.  Returns true on unrecoverable failure.
  bool advance(double t_start, double t_end);

 protected:
  // Hooks called after each inner step — override for custom behavior.
  // Default onSuccess_() calls vis, observations, and checkpoint.
  // Default onFail_() calls failed visualization.
  virtual void onSuccess_(double t_old, double t_new);
  virtual void onFail_(double t_old, double t_new);

 private:
  // True when the loop should stop (besides reaching t_end).
  bool extraDoneCheck_(double t, int cycle) const;

  // vis/obs/checkpoint helpers
  bool visualize_(bool force = false);
  void observe_();
  bool checkpoint_(bool force = false);

  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<Amanzi::VerboseObject> vo_;
  Teuchos::RCP<Amanzi::State> S_;
  Teuchos::RCP<Amanzi::PK> pk_;
  Teuchos::RCP<Amanzi::Utils::TimeStepManager> tsm_;
  Amanzi::Tag tag_current_, tag_next_;
  Teuchos::RCP<Teuchos::Time> wallclock_timer_;

  double max_dt_, min_dt_;
  int cycle1_;
  double duration_;
  bool subcycled_ts_;

  std::vector<Teuchos::RCP<Amanzi::Visualization>> visualization_;
  std::vector<Teuchos::RCP<Amanzi::Visualization>> failed_visualization_;
  Teuchos::RCP<Amanzi::Checkpoint> checkpoint_obj_;
  std::vector<Teuchos::RCP<Amanzi::UnstructuredObservations>> observations_;
};

} // namespace ATS
