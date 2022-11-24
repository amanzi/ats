#include "liquid_gas_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, LiquidGasWaterContentEvaluator>
  LiquidGasWaterContentEvaluator::reg_("liquid+gas water content");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
