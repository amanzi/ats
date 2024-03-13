/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "albedo_twocomponent_model.hh"
#include "albedo_threecomponent_model.hh"
#include "area_fractions_twocomponent_model.hh"
#include "area_fractions_threecomponent_model.hh"
// #include "area_fractions_threecomponent_microtopography_evaluator.hh"
#include "canopy_drainage_model.hh"
#include "evaporation_downregulation_soil_model.hh"
#include "incident_shortwave_radiation_evaluator.hh"
#include "incoming_longwave_radiation_model.hh"
#include "interception_fraction_model.hh"
#include "pet_priestley_taylor_evaluator.hh"
// #include "plant_wilting_factor_evaluator.hh"
#include "rooting_depth_fraction_evaluator.hh"
#include "snow_meltrate_model.hh"
// #include "transpiration_distribution_evaluator.hh"
#include "transpiration_distribution_relperm_evaluator.hh"
#include "radiation_balance_evaluator.hh"
// #include "canopy_radiation_evaluator.hh"
// #include "seb_twocomponent_evaluator.hh"
// #include "seb_threecomponent_evaluator.hh"

namespace Amanzi {

template<>
REGISTER(SurfaceBalance::Relations::AlbedoTwoComponentEvaluator);
template<>
REGISTER(SurfaceBalance::Relations::AlbedoThreeComponentEvaluator);

template<>
REGISTER(SurfaceBalance::Relations::AreaFractionsTwoComponentEvaluator);
template<>
REGISTER(SurfaceBalance::Relations::AreaFractionsThreeComponentEvaluator);

REGISTER(SurfaceBalance::Relations::IncidentShortwaveRadiationEvaluator);

REGISTER_MODEL(SurfaceBalance::Relations::IncomingLongwaveRadiationModel);
REGISTER_MODEL(SurfaceBalance::Relations::InterceptionFractionModel);
REGISTER_MODEL(SurfaceBalance::Relations::CanopyDrainageModel);

REGISTER(SurfaceBalance::Relations::PETPriestleyTaylorEvaluator);
REGISTER_MODEL(SurfaceBalance::Relations::EvaporationDownregulationSoilModel);

// Utils::RegisteredFactory<Evaluator, PlantWiltingFactorEvaluator>
//   PlantWiltingFactorEvaluator::reg_("plant wilting factor");

REGISTER(SurfaceBalance::Relations::RootingDepthFractionEvaluator);

// Utils::RegisteredFactory<Evaluator, TranspirationDistributionEvaluator>
//   TranspirationDistributionEvaluator::reg_("transpiration distribution via rooting depth");
REGISTER(SurfaceBalance::Relations::TranspirationDistributionRelPermEvaluator);

REGISTER_BY_MATERIAL(SurfaceBalance::Relations::SnowMeltRateModel);

REGISTER(SurfaceBalance::Relations::RadiationBalanceEvaluator);

// Utils::RegisteredFactory<Evaluator, CanopyRadiationEvaluator>
//   CanopyRadiationEvaluator::reg_("canopy radiation balance from above");

// Utils::RegisteredFactory<Evaluator, SEBTwoComponentEvaluator>
//   SEBTwoComponentEvaluator::reg_("surface energy balance, two components");

// Utils::RegisteredFactory<Evaluator, SEBThreeComponentEvaluator>
//   SEBThreeComponentEvaluator::reg_("surface energy balance, three components");

} // namespace Amanzi
