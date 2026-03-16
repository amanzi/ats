/*
  ATS-EcoSIM, Code Adapted for use from Alquimia

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson

  This is a point of contact for the chemistry engine exposed by Alquimia
  to the rest of Amanzi--it provides the ability to enforce geochemical
  conditions and to integrate reactions given a chemical configuration.
*/

#ifndef BGC_ENGINE_HH_
#define BGC_ENGINE_HH_

#include <string>
#include <vector>
#include <map>

#include "BGC_memory.hh"
#include "BGC_containers.hh"

#include "VerboseObject.hh"

namespace Amanzi {
namespace EcoSIM {

class BGCEngine {
 public:

  // Constructs a chemistry engine using the given engine (backend) name and input file.
  BGCEngine(const std::string& engineName, const std::string& inputFile);

  // Destructor.
  ~BGCEngine();

  // Returns the name of the backend that does the chemistry.
  const std::string& Name() const;

  // Returns true if the chemistry engine is thread-safe, false if not.
  bool IsThreadSafe() const;

  // Returns a reference to a "sizes" object that can be queried to find the sizes of the various
  // arrays representing the geochemical state within the engine.
  const BGCSizes& Sizes() const;

  // Initializes the data structures that hold the chemical state information.
  void InitState(BGCProperties& properties,
                 BGCState& state,
                 BGCAuxiliaryData& aux_data,
                 int ncells_per_col_,
                 int num_components,
                 int num_columns);

  // Frees the data structures that hold the chemical state information.
  void FreeState(BGCProperties& properties,
                 BGCState& state,
                 BGCAuxiliaryData& aux_data);

  void DataTest();

  bool Setup(BGCProperties& properties,
               BGCState& state,
               BGCSizes& sizes,
               int num_iterations,
               int num_columns,
               int ncells_per_col_);

  bool Advance(const double delta_time,
               BGCProperties& properties,
               BGCState& state,
               BGCSizes& sizes,
               int num_iterations,
               int num_columns);

  void CopyBGCState(const BGCState* const source,
                         BGCState* destination);
  void CopyBGCProperties(const BGCProperties* const source,
                              BGCProperties* destination);

 private:

  // bgc data structures.
  bool bgc_initialized_;
  void* engine_state_;
  BGCSizes sizes_;
  BGCInterface bgc_;

  Teuchos::RCP<VerboseObject> vo_;
  // Back-end engine name and input file.
  std::string bgc_engine_name_;
  std::string bgc_engine_inputfile_;

  // forbidden.
  BGCEngine();
  BGCEngine(const BGCEngine&);
  BGCEngine& operator=(const BGCEngine&);

};

} // namespace
} // namespace

#endif
