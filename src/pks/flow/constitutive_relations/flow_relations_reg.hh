#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "richards_water_content_model.hh"
#include "compressible_porosity_linear_model.hh"
#include "capillary_pressure_liquid_atm_model.hh"

#include "wrm_van_genuchten.hh"
#include "wrm_model.hh"
#include "relative_permeability_model.hh"

namespace Amanzi {

REGISTER(Flow::Relations::RichardsWaterContentModel);
REGISTER(Flow::Relations::CompressiblePorosityLinearModel);
REGISTER(Flow::Relations::CapillaryPressureLiquidAtmModel);

REGISTER(Flow::Relations::WRMVanGenuchtenModel);
REGISTER(Flow::Relations::RelativePermeabilityVanGenuchtenModel);


} // namespace Amanzi
