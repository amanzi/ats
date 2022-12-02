/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

#include "thermal_conductivity_twophase_wetdry.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityTwoPhase, ThermalConductivityTwoPhaseWetDry>
  ThermalConductivityTwoPhaseWetDry::factory_("two-phase wet/dry");

} // namespace Energy
} // namespace Amanzi
