#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "compressible_porosity_linear_model.hh"

namespace Amanzi {

REGISTER_MODEL(Flow::Relations::CompressiblePorosityLinearModel);
REGISTER_BY_MATERIAL(Flow::Relations::CompressiblePorosityLinearModel);

} // namespace Amanzi
