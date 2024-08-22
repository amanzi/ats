#include "density_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<Evaluator,DensityEvaluator> DensityEvaluator::factory_("density");

} //namespace
} //namespace
