#include "evaporative_flux_relaxation_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

Utils::RegisteredFactory<Evaluator, EvaporativeFluxRelaxationEvaluator>
  EvaporativeFluxRelaxationEvaluator::reg_("evaporative flux relaxation");

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
