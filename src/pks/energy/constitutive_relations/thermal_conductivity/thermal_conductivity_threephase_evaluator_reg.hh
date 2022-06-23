/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------

   ATS

   Three-phase TC evaluator
   ------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_evaluator.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<Evaluator,ThermalConductivityThreePhaseEvaluator> 
ThermalConductivityThreePhaseEvaluator::factory_("three-phase thermal conductivity");

} // namespace
} // namespace
