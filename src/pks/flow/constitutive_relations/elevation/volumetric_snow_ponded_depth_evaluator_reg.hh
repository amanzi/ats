#include "volumetric_snow_ponded_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,VolumetricSnowPondedDepthEvaluator> VolumetricSnowPondedDepthEvaluator::reg_("volumetric ponded and snow depths");

} //namespace
} //namespace
