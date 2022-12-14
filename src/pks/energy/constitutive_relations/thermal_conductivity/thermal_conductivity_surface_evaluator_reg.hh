/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------

   ATS

   Surface TC evaluator
   ------------------------------------------------------------------------- */

#include "thermal_conductivity_surface_evaluator.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<Evaluator, ThermalConductivitySurfaceEvaluator>
  ThermalConductivitySurfaceEvaluator::factory_("surface thermal conductivity");

} // namespace Energy
} // namespace Amanzi
