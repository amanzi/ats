/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/* -------------------------------------------------------------------------

   ATS

   Surface TC evaluator
   ------------------------------------------------------------------------- */

#include "thermal_conductivity_surface_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

// registry of method
Utils::RegisteredFactory<Evaluator, ThermalConductivitySurfaceEvaluator>
  ThermalConductivitySurfaceEvaluator::factory_("surface thermal conductivity");

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
