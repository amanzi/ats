/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_volume_averaged.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityThreePhase, ThermalConductivityThreePhaseVolumeAveraged>
  ThermalConductivityThreePhaseVolumeAveraged::factory_("three-phase volume averaged");


} // namespace Energy
} // namespace Amanzi
