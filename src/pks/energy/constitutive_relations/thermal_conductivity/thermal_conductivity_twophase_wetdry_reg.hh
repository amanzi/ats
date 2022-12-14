/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon


Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.

See ATS process model documentation's permafrost model for details.
------------------------------------------------------------------------- */


#include "thermal_conductivity_twophase_wetdry.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityTwoPhase, ThermalConductivityTwoPhaseWetDry>
  ThermalConductivityTwoPhaseWetDry::factory_("two-phase wet/dry");

} // namespace Energy
} // namespace Amanzi
