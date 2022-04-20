/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


#include "albedo_twocomponent_evaluator.hh"
#include "albedo_threecomponent_evaluator.hh"
#include "area_fractions_twocomponent_evaluator.hh"
#include "area_fractions_threecomponent_evaluator.hh"
#include "area_fractions_threecomponent_microtopography_evaluator.hh"
#include "drainage_evaluator.hh"
#include "evaporation_downregulation_evaluator.hh"
#include "incident_shortwave_radiation_evaluator.hh"
#include "longwave_evaluator.hh"
#include "interception_fraction_evaluator.hh"
#include "pet_priestley_taylor_evaluator.hh"
#include "plant_wilting_factor_evaluator.hh"
#include "rooting_depth_fraction_evaluator.hh"
#include "snow_meltrate_evaluator.hh"
#include "transpiration_distribution_evaluator.hh"
#include "radiation_balance_evaluator.hh"
#include "seb_twocomponent_evaluator.hh"
#include "seb_threecomponent_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

Utils::RegisteredFactory<FieldEvaluator,AlbedoTwoComponentEvaluator> AlbedoTwoComponentEvaluator::reg_("subgrid albedos, two components");

Utils::RegisteredFactory<FieldEvaluator,AlbedoThreeComponentEvaluator> AlbedoThreeComponentEvaluator::reg_("subgrid albedos, three components");

Utils::RegisteredFactory<FieldEvaluator,AreaFractionsTwoComponentEvaluator> AreaFractionsTwoComponentEvaluator::reg_("area fractions, two components");

Utils::RegisteredFactory<FieldEvaluator,AreaFractionsThreeComponentEvaluator> AreaFractionsThreeComponentEvaluator::reg_("area fractions, three components");

Utils::RegisteredFactory<FieldEvaluator,AreaFractionsThreeComponentMicrotopographyEvaluator> AreaFractionsThreeComponentMicrotopographyEvaluator::reg_("area fractions, three components with microtopography");

Utils::RegisteredFactory<FieldEvaluator,IncidentShortwaveRadiationEvaluator> IncidentShortwaveRadiationEvaluator::reg_("incident shortwave radiation");

Utils::RegisteredFactory<FieldEvaluator,LongwaveEvaluator> LongwaveEvaluator::reg_("incoming longwave radiation");

Utils::RegisteredFactory<FieldEvaluator,InterceptionFractionEvaluator> InterceptionFractionEvaluator::reg_("interception fraction");

Utils::RegisteredFactory<FieldEvaluator,DrainageEvaluator> DrainageEvaluator::reg_("canopy drainage");

Utils::RegisteredFactory<FieldEvaluator,PETPriestleyTaylorEvaluator> PETPriestleyTaylorEvaluator::reg_("potential evapotranspiration, Priestley-Taylor");

Utils::RegisteredFactory<FieldEvaluator,EvaporationDownregulationEvaluator> EvaporationDownregulationEvaluator::reg_("evaporation downregulation via soil resistance");

Utils::RegisteredFactory<FieldEvaluator,PlantWiltingFactorEvaluator> PlantWiltingFactorEvaluator::reg_("plant wilting factor");

Utils::RegisteredFactory<FieldEvaluator,RootingDepthFractionEvaluator> RootingDepthFractionEvaluator::reg_("rooting depth fraction");

Utils::RegisteredFactory<FieldEvaluator,TranspirationDistributionEvaluator> TranspirationDistributionEvaluator::reg_("transpiration distribution via rooting depth");

Utils::RegisteredFactory<FieldEvaluator,SnowMeltRateEvaluator> SnowMeltRateEvaluator::reg_("snow melt rate");

Utils::RegisteredFactory<FieldEvaluator,RadiationBalanceEvaluator> RadiationBalanceEvaluator::reg_("radiation balance, surface and canopy");

Utils::RegisteredFactory<FieldEvaluator,SEBTwoComponentEvaluator> SEBTwoComponentEvaluator::reg_("surface energy balance, two components");

Utils::RegisteredFactory<FieldEvaluator,SEBThreeComponentEvaluator> SEBThreeComponentEvaluator::reg_("surface energy balance, three components");

} // namespace
} // namespace
} // namespace
