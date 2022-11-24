#include "liquid_ice_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, LiquidIceWaterContentEvaluator>
  LiquidIceWaterContentEvaluator::reg_("liquid+ice water content");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
