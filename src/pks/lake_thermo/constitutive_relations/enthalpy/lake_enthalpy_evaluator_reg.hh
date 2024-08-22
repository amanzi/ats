#include "lake_enthalpy_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<Evaluator,LakeEnthalpyEvaluator> LakeEnthalpyEvaluator::factory_("lake enthalpy");

} //namespace
} //namespace
