#include "bioturbation_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

// registry of method
Utils::RegisteredFactory<Evaluator, BioturbationEvaluator>
  BioturbationEvaluator::fac_("bioturbation evaluator");

} // namespace BGCRelations
} // namespace BGC
} // namespace Amanzi
