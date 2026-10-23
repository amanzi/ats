/*
  ATS-EcoSIM, Code Adapted for use from Alquimia

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson

  Point of contact between the EcoSIM PK and EcoSIM: allocates the exchange
  containers (EcoContainers.hh) and calls EcoSIM's C entry points
  (ecosim_interface.h).
*/

#ifndef ECO_ENGINE_HH_
#define ECO_ENGINE_HH_

#include <string>
#include <vector>

#include "EcoMemory.hh"
#include "EcoContainers.hh"

namespace Amanzi {
namespace EcoSIM {

// One entry of the EcoSIM internal state layout (see ATSStateRegistryMod.F90)
struct EcoInternalStateEntry {
  std::string ats_name;    // ATS field basename, key "surface-<ats_name>"
  std::string ecosim_name; // EcoSIM variable holding it
  std::string units;       // as annotated in EcoSIM
  std::string description;
  int num_components;      // per column
  int role;                // kEcoRolePrivate or kEcoRoleOutput
};

class EcoEngine {
 public:
  EcoEngine(const std::string& engine_name, const std::string& input_file);
  ~EcoEngine();

  // Allocates/frees the arrays of the environment and feedback containers.
  void InitState(EcoEnvironment& environment,
                 EcoFeedback& feedback,
                 int ncells_per_col_,
                 int num_components,
                 int num_columns,
                 int num_pfts);
  void FreeState(EcoEnvironment& environment, EcoFeedback& feedback);

  // Throws if the container sizes differ between this build and EcoSIM's
  // bind(C) types (a mismatch means the C header and EcoContainers.F90 differ).
  void CheckContainerSizes() const;

  void DataTest();

  bool Setup(EcoConfig& config,
             EcoEnvironment& environment,
             EcoFeedback& feedback,
             EcoInternalState& internal_state,
             EcoSizes& sizes,
             int num_iterations,
             int num_columns,
             int ncells_per_col_);

  bool Advance(const double delta_time,
               EcoEnvironment& environment,
               EcoFeedback& feedback,
               EcoInternalState& internal_state,
               EcoSizes& sizes,
               int num_iterations,
               int num_columns);

  // Layout of the EcoSIM internal state. Static: valid before Setup.
  int InternalStateLayoutVersion() const;
  std::vector<EcoInternalStateEntry> InternalStateLayout(const EcoSizes& sizes) const;

  // Allocates/frees the container for the internal state of num_columns columns.
  void InitInternalState(EcoInternalState& internal_state,
                         const std::vector<EcoInternalStateEntry>& layout,
                         int num_columns);
  void FreeInternalState(EcoInternalState& internal_state);

 private:
  EcoInterface eco_;
  std::string engine_name_;
  std::string engine_inputfile_;

  // forbidden.
  EcoEngine();
  EcoEngine(const EcoEngine&);
  EcoEngine& operator=(const EcoEngine&);
};

} // namespace EcoSIM
} // namespace Amanzi

#endif
