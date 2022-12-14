#include "overland_conductivity_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, OverlandConductivityEvaluator>
  OverlandConductivityEvaluator::factory_("overland conductivity");

} // namespace Flow
} // namespace Amanzi
