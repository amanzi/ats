#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "richards_water_content_model.hh"

namespace Amanzi {

REGISTER(Flow::Relations::RichardsWaterContentModel);

} // namespace Amanzi
