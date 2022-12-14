#include "icy_height_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, IcyHeightEvaluator>
  IcyHeightEvaluator::factory_("icy ponded depth");

} // namespace Flow
} // namespace Amanzi
