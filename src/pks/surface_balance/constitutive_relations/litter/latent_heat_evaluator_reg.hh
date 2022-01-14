#include "latent_heat_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

Utils::RegisteredFactory<Evaluator,LatentHeatEvaluator> LatentHeatEvaluator::reg_("latent heat from evaporative flux");

} //namespace
} //namespace
} //namespace
