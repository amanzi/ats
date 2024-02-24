#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "capillary_pressure_liquid_atm_model.hh"

#include "wrm_van_genuchten.hh"
#include "wrm_model.hh"
#include "relative_permeability_model.hh"
#include "relative_permeability_evaluator.hh"

namespace Amanzi {

REGISTER_MODEL(Flow::Relations::CapillaryPressureLiquidAtmModel);

REGISTER_MODEL(Flow::Relations::WRMVanGenuchtenModel);
REGISTER_BY_MATERIAL(Flow::Relations::WRMVanGenuchtenModel);


template<>
REGISTER(Flow::Relations::RelativePermeabilityEvaluator<EvaluatorModelCV<RelativePermeabilityVanGenuchtenModel>>);

template<>
REGISTER(Flow::Relations::RelativePermeabilityEvaluator<EvaluatorModelCVByMaterial<RelativePermeabilityVanGenuchtenModel>>);


} // namespace Amanzi
