#include "pool_decomposition_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

// registry of method
Utils::RegisteredFactory<Evaluator, PoolDecompositionEvaluator>
  PoolDecompositionEvaluator::fac_("pool decomposition evaluator");

} // namespace BGCRelations
} // namespace BGC
} // namespace Amanzi
