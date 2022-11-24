#include "distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, DistributedTilesRateEvaluator>
  DistributedTilesRateEvaluator::reg_("distributed tiles");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
