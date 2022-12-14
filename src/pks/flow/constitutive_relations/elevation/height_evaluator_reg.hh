#include "height_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, HeightEvaluator> HeightEvaluator::factory_("ponded depth");

} // namespace Flow
} // namespace Amanzi
