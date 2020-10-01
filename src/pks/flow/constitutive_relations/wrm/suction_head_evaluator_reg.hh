#include "suction_head_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SuctionHeadEvaluator> SuctionHeadEvaluator::factory_("WRM suction head");

} //namespace
} //namespace
