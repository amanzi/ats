/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "energy_three_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<ThreePhase> ThreePhase::reg_("three-phase energy");

} // namespace Energy
} // namespace Amanzi
