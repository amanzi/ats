#include "suction_head_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, SuctionHeadEvaluator>
  SuctionHeadEvaluator::factory_("WRM suction head");

} // namespace Flow
} // namespace Amanzi
