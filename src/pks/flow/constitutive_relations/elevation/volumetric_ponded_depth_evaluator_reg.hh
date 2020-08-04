#include "volumetric_ponded_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,VolumetricPondedDepthEvaluator> VolumetricPondedDepthEvaluator::reg_("volumetric ponded depth");

} //namespace
} //namespace
