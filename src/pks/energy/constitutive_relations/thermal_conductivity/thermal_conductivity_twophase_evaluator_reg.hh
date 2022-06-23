/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------

   ATS

   Two-phase TC evaluator
   ------------------------------------------------------------------------- */

#include "thermal_conductivity_twophase_evaluator.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<Evaluator,ThermalConductivityTwoPhaseEvaluator> 
ThermalConductivityTwoPhaseEvaluator::factory_("two-phase thermal conductivity");

} // namespace
} // namespace
