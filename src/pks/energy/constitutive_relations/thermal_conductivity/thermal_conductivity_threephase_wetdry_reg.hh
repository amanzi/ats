/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_wetdry.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityThreePhase, ThermalConductivityThreePhaseWetDry>
  ThermalConductivityThreePhaseWetDry::factory_("three-phase wet/dry");

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
