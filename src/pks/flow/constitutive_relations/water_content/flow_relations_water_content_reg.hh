#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "richards_water_content_model.hh"
#include "overland_pressure_water_content_model.hh"
#include "ponded_depth_model.hh"

namespace Amanzi {

REGISTER_MODEL(Flow::Relations::RichardsWaterContentModel);
REGISTER_MODEL(Flow::Relations::OverlandPressureWaterContentModel);
REGISTER_MODEL(Flow::Relations::PondedDepthModel);

} // namespace Amanzi
