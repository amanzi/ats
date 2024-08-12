#include "lake_heat_capacity_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,LakeHeatCapacityEvaluator> LakeHeatCapacityEvaluator::factory_("lake heat capacity");

} //namespace
} //namespace
