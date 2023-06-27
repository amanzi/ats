#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "capillary_pressure_liquid_atm_model.hh"

#include "wrm_van_genuchten.hh"
#include "wrm_model.hh"
#include "relative_permeability_model.hh"

namespace Amanzi {

REGISTER(Flow::Relations::CapillaryPressureLiquidAtmModel);

REGISTER(Flow::Relations::WRMVanGenuchtenModel);
REGISTER_BY_MATERIAL(Flow::Relations::WRMVanGenuchtenModel);

REGISTER(Flow::Relations::RelativePermeabilityVanGenuchtenModel);
REGISTER_BY_MATERIAL(Flow::Relations::RelativePermeabilityVanGenuchtenModel);


} // namespace Amanzi
