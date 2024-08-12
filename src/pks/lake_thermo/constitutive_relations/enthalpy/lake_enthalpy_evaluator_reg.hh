#include "lake_enthalpy_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,LakeEnthalpyEvaluator> LakeEnthalpyEvaluator::factory_("lake enthalpy");

} //namespace
} //namespace
