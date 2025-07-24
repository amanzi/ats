/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A set of helper functions for doing common things in PKs.
#pragma once

#include <string>
#include "Key.hh"

namespace Amanzi {

class State;

// -----------------------------------------------------------------------------
// Helper functions for working with Amanzi's Chemistry PK
// -----------------------------------------------------------------------------
void convertConcentrationToMolFrac(State& S,
                                   const KeyTag& tcc,
                                   const KeyTag& mol_frac,
                                   const KeyTag& mol_dens,
                                   const std::string& passwd);


void convertMolFracToConcentration(State& S,
                                   const KeyTag& mol_frac,
                                   const KeyTag& tcc,
                                   const KeyTag& mol_dens,
                                   const std::string& passwd);

} // namespace Amanzi
