#include "snow_skin_potential_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, SnowSkinPotentialEvaluator>
  SnowSkinPotentialEvaluator::factory_("snow skin potential");

} // namespace Flow
} // namespace Amanzi
