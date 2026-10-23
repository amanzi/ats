/*
  Basic architecture based on the Alquima interfece adapted for
  use in the ATSEcoSIM PK

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson
           Sergi Molins <smolins@lbl.gov>

  Calls into EcoSIM through the exchange containers.
*/

#include <cstring>
#include <sstream>
#include "EcoEngine.hh"
#include "ecosim_interface.h"
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {
namespace EcoSIM {

EcoEngine::EcoEngine(const std::string& engine_name, const std::string& input_file)
  : engine_name_(engine_name), engine_inputfile_(input_file)
{
  CreateEcoInterface(engine_name_.c_str(), &eco_);
}

EcoEngine::~EcoEngine()
{
  eco_.Shutdown();
}

void EcoEngine::InitState(EcoEnvironment& environment,
                          EcoFeedback& feedback,
                          int ncells_per_col_,
                          int num_components,
                          int num_columns,
                          int num_pfts)
{
  AllocateEcoEnvironment(&environment, ncells_per_col_, num_components, num_columns);
  AllocateEcoFeedback(&feedback, ncells_per_col_, num_columns, num_pfts);
}

void EcoEngine::FreeState(EcoEnvironment& environment, EcoFeedback& feedback)
{
  FreeEcoEnvironment(&environment);
  FreeEcoFeedback(&feedback);
}

void EcoEngine::CheckContainerSizes() const
{
  // same order as ecosim_container_sizes() in EcoSIM's ecosim_wrappers.F90
  const char* names[kEcoNumContainerTypes] = { "EcoVectorDouble", "EcoMatrixDouble",
                                               "EcoTensorDouble", "EcoSizes",
                                               "EcoConfig",       "EcoEnvironment",
                                               "EcoFeedback",     "EcoInternalState" };
  const size_t ats_sizes[kEcoNumContainerTypes] = {
    sizeof(EcoVectorDouble), sizeof(EcoMatrixDouble), sizeof(EcoTensorDouble),
    sizeof(EcoSizes),        sizeof(EcoConfig),       sizeof(EcoEnvironment),
    sizeof(EcoFeedback),     sizeof(EcoInternalState)
  };
  size_t ecosim_sizes[kEcoNumContainerTypes] = { 0 };
  int n = kEcoNumContainerTypes;
  ecosim_container_sizes(&n, ecosim_sizes);

  std::stringstream bad;
  if (n != kEcoNumContainerTypes)
    bad << "  EcoSIM reports " << n << " container types, ATS expects "
        << kEcoNumContainerTypes << "\n";
  for (int i = 0; i != kEcoNumContainerTypes && i != n; ++i) {
    if (ats_sizes[i] != ecosim_sizes[i])
      bad << "  " << names[i] << ": " << ats_sizes[i] << " bytes in ATS, " << ecosim_sizes[i]
          << " bytes in EcoSIM\n";
  }
  if (!bad.str().empty()) {
    Errors::Message msg;
    msg << "EcoSIM: the exchange containers differ between ATS (data/EcoContainers.hh) and "
        << "EcoSIM (ATSUtils/EcoContainers.F90):\n"
        << bad.str();
    Exceptions::amanzi_throw(msg);
  }
}

// The internal state layout is defined by EcoSIM, so it is queried directly
// rather than through the EcoInterface function table.
int EcoEngine::InternalStateLayoutVersion() const
{
  return ecosim_internal_state_layout_version();
}

std::vector<EcoInternalStateEntry>
EcoEngine::InternalStateLayout(const EcoSizes& sizes) const
{
  std::vector<EcoInternalStateEntry> layout;
  int num_entries = ecosim_internal_state_num_entries(&sizes);
  for (int i = 0; i != num_entries; ++i) {
    char ats_name[kEcoStateNameLength], ecosim_name[kEcoStateNameLength];
    char units[kEcoStateUnitsLength], description[kEcoStateDescriptionLength];
    EcoInternalStateEntry entry;
    ecosim_internal_state_entry(
      i, &sizes, ats_name, ecosim_name, units, description, &entry.num_components, &entry.role);
    entry.ats_name = ats_name;
    entry.ecosim_name = ecosim_name;
    entry.units = units;
    entry.description = description;
    layout.push_back(entry);
  }
  return layout;
}

void EcoEngine::InitInternalState(EcoInternalState& internal_state,
                                  const std::vector<EcoInternalStateEntry>& layout,
                                  int num_columns)
{
  int values_per_column = 0;
  for (const auto& entry : layout) values_per_column += entry.num_components;
  AllocateEcoInternalState(&internal_state, InternalStateLayoutVersion(), layout.size(),
                           num_columns, values_per_column);
}

void EcoEngine::FreeInternalState(EcoInternalState& internal_state)
{
  FreeEcoInternalState(&internal_state);
}

void EcoEngine::DataTest()
{
  eco_.DataTest();
}

bool EcoEngine::Setup(EcoConfig& config,
                      EcoEnvironment& environment,
                      EcoFeedback& feedback,
                      EcoInternalState& internal_state,
                      EcoSizes& sizes,
                      int num_iterations,
                      int num_columns,
                      int ncells_per_col_)
{
  eco_.Setup(&config, &environment, &feedback, &internal_state, &sizes, num_iterations,
             num_columns, ncells_per_col_);
  return true;
}

bool EcoEngine::Advance(const double delta_time,
                        EcoEnvironment& environment,
                        EcoFeedback& feedback,
                        EcoInternalState& internal_state,
                        EcoSizes& sizes,
                        int num_iterations,
                        int num_columns)
{
  eco_.Advance(delta_time, &environment, &feedback, &internal_state, &sizes, num_iterations,
               num_columns);
  return true;
}

} // namespace EcoSIM
} // namespace Amanzi
