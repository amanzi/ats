/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include "thermal_conductivity_twophase_peterslidard.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityTwoPhase, ThermalConductivityTwoPhasePetersLidard>
  ThermalConductivityTwoPhasePetersLidard::factory_("two-phase Peters-Lidard");

} // namespace Energy
} // namespace Amanzi
