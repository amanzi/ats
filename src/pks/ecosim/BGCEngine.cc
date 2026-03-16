/*
  Basic architecture based on the Alquima interfece adapted for
  use in the ATSEcoSIM PK

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson
           Sergi Molins <smolins@lbl.gov>

  This implements the Alquimia chemistry engine.
*/

#include <iostream>
#include <cstring>
#include <cstdio>
#include <assert.h>
#include "BGCEngine.hh"
#include "errors.hh"
#include "exceptions.hh"

// Support for manipulating floating point exception handling.
#ifdef _GNU_SOURCE
#define AMANZI_USE_FENV
#include <fenv.h>
#endif

namespace Amanzi {
namespace EcoSIM {

BGCEngine::BGCEngine(const std::string& engineName,
                                 const std::string& inputFile) :
  bgc_engine_name_(engineName),
  bgc_engine_inputfile_(inputFile)
{
  Errors::Message msg;

  CreateBGCInterface(bgc_engine_name_.c_str(),
                    &bgc_);

}

BGCEngine::~BGCEngine()
{
  bgc_.Shutdown();

  //Did I forget to implement this?
  //FreeBGCProperties(&props);
  //FreeBGCState(&state);
  //FreeBGCAuxiliaryData(&aux_data);
  //FreeAlquimiaEngineStatus(&chem_status_);
}

const BGCSizes&
BGCEngine::Sizes() const
{
  return sizes_;
}

void BGCEngine::InitState(BGCProperties& properties,
                                BGCState& state,
                                BGCAuxiliaryData& aux_data,
                                int ncells_per_col_,
                                int num_components,
                                int num_columns)
{
  AllocateBGCProperties(&sizes_, &properties, ncells_per_col_, num_columns);
  AllocateBGCState(&sizes_, &state, ncells_per_col_, num_components, num_columns);
}

void BGCEngine::FreeState(BGCProperties& properties,
                                BGCState& state,
                                BGCAuxiliaryData& aux_data)
{
  FreeBGCProperties(&properties);
  FreeBGCState(&state);
}

void BGCEngine::DataTest() {

  bgc_.DataTest();
}

bool BGCEngine::Setup(BGCProperties& properties,
                              BGCState& state,
                              BGCSizes& sizes_,
                              int num_iterations,
                              int num_columns,
                              int ncells_per_col_)
{
  bgc_.Setup(&properties,
                &state,
                &sizes_,
                num_iterations,
                num_columns,
                ncells_per_col_);

}

bool BGCEngine::Advance(const double delta_time,
                              BGCProperties& properties,
                              BGCState& state,
                              BGCSizes& sizes_,
                              int num_iterations,
                              int num_columns)
{
  bgc_.Advance(delta_time,
                &properties,
                &state,
                &sizes_,
                num_iterations,
                num_columns);

}

} // namespace
} // namespace
