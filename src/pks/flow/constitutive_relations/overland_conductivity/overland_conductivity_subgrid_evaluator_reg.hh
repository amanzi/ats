#include "overland_conductivity_subgrid_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,OverlandConductivitySubgridEvaluator> OverlandConductivitySubgridEvaluator::factory_("overland conductivity subgrid");

}
}
