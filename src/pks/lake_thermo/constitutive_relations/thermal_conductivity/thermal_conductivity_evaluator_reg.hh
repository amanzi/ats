#include "thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ThermalConductivityEvaluator> ThermalConductivityEvaluator::factory_("lake thermal conductivity");

} //namespace
} //namespace
