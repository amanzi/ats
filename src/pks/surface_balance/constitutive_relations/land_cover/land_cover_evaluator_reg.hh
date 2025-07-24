/*
  Copyright 2010-202x held jointly by participating institutions.
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
#include "transpiration_distribution_relperm_evaluator.hh"
#include "radiation_balance_evaluator.hh"
#include "canopy_radiation_evaluator.hh"
#include "seb_twocomponent_evaluator.hh"
#include "seb_threecomponent_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

Utils::RegisteredFactory<Evaluator, AlbedoTwoComponentEvaluator> AlbedoTwoComponentEvaluator::reg_(
  "subgrid albedos, two components");

Utils::RegisteredFactory<Evaluator, AlbedoThreeComponentEvaluator>
  AlbedoThreeComponentEvaluator::reg_("subgrid albedos, three components");

Utils::RegisteredFactory<Evaluator, AreaFractionsTwoComponentEvaluator>
  AreaFractionsTwoComponentEvaluator::reg_("area fractions, two components");

Utils::RegisteredFactory<Evaluator, AreaFractionsThreeComponentEvaluator>
  AreaFractionsThreeComponentEvaluator::reg_("area fractions, three components");

Utils::RegisteredFactory<Evaluator, AreaFractionsThreeComponentMicrotopographyEvaluator>
  AreaFractionsThreeComponentMicrotopographyEvaluator::reg_(
    "area fractions, three components with microtopography");

Utils::RegisteredFactory<Evaluator, IncidentShortwaveRadiationEvaluator>
  IncidentShortwaveRadiationEvaluator::reg_("incident shortwave radiation");

Utils::RegisteredFactory<Evaluator, LongwaveEvaluator> LongwaveEvaluator::reg_(
  "incoming longwave radiation");

Utils::RegisteredFactory<Evaluator, InterceptionFractionEvaluator>
  InterceptionFractionEvaluator::reg_("interception fraction");

Utils::RegisteredFactory<Evaluator, DrainageEvaluator> DrainageEvaluator::reg_("canopy drainage");

Utils::RegisteredFactory<Evaluator, PETPriestleyTaylorEvaluator> PETPriestleyTaylorEvaluator::reg_(
  "potential evapotranspiration, Priestley-Taylor");

Utils::RegisteredFactory<Evaluator, EvaporationDownregulationEvaluator>
  EvaporationDownregulationEvaluator::reg_("evaporation downregulation, soil resistance");

Utils::RegisteredFactory<Evaluator, PlantWiltingFactorEvaluator> PlantWiltingFactorEvaluator::reg_(
  "plant wilting factor");

Utils::RegisteredFactory<Evaluator, RootingDepthFractionEvaluator>
  RootingDepthFractionEvaluator::reg_("root fraction");

Utils::RegisteredFactory<Evaluator, TranspirationDistributionEvaluator>
  TranspirationDistributionEvaluator::reg_("transpiration distribution, rooting depth");
Utils::RegisteredFactory<Evaluator, TranspirationDistributionRelPermEvaluator>
  TranspirationDistributionRelPermEvaluator::reg_(
    "transpiration distribution, relative permeability");

Utils::RegisteredFactory<Evaluator, SnowMeltRateEvaluator> SnowMeltRateEvaluator::reg_(
  "snow melt rate");

Utils::RegisteredFactory<Evaluator, RadiationBalanceEvaluator> RadiationBalanceEvaluator::reg_(
  "radiation balance, surface and canopy");

Utils::RegisteredFactory<Evaluator, CanopyRadiationEvaluator> CanopyRadiationEvaluator::reg_(
  "canopy radiation balance from above");

Utils::RegisteredFactory<Evaluator, SEBTwoComponentEvaluator> SEBTwoComponentEvaluator::reg_(
  "surface energy balance, two components");

Utils::RegisteredFactory<Evaluator, SEBThreeComponentEvaluator> SEBThreeComponentEvaluator::reg_(
  "surface energy balance, three components");

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
