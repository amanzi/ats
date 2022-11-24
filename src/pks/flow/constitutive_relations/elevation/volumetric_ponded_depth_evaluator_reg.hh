#include "volumetric_ponded_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, VolumetricPondedDepthEvaluator>
  VolumetricPondedDepthEvaluator::reg_("volumetric ponded depth");

} // namespace Flow
} // namespace Amanzi
