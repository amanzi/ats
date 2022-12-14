/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_sutra_hacked.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityThreePhase, ThermalConductivityThreePhaseSutraHacked>
  ThermalConductivityThreePhaseSutraHacked::factory_("three-phase sutra hacked");


} // namespace Energy
} // namespace Amanzi
