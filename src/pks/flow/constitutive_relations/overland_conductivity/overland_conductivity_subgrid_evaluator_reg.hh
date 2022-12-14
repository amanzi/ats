#include "overland_conductivity_subgrid_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, OverlandConductivitySubgridEvaluator>
  OverlandConductivitySubgridEvaluator::factory_("overland conductivity subgrid");

} // namespace Flow
} // namespace Amanzi
