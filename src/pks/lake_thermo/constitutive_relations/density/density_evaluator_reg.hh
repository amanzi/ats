#include "density_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,DensityEvaluator> DensityEvaluator::factory_("density");

} //namespace
} //namespace
