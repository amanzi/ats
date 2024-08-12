#include "lake_energy_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,LakeEnergyEvaluator> LakeEnergyEvaluator::factory_("lake energy");

} //namespace
} //namespace
