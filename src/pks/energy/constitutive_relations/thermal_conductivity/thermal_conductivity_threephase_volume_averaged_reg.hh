/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "thermal_conductivity_threephase_volume_averaged.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityThreePhase, ThermalConductivityThreePhaseVolumeAveraged>
  ThermalConductivityThreePhaseVolumeAveraged::factory_("three-phase volume averaged");


} // namespace Energy
} // namespace Amanzi
