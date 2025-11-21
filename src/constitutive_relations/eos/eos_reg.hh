/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "eos_evaluator.hh"
#include "eos_constant.hh"
#include "eos_ice.hh"
#include "eos_ideal_gas.hh"
#include "eos_linear.hh"
#include "eos_sw.hh"
#include "eos_vapor_in_gas.hh"
#include "eos_water.hh"

#include "molar_fraction_gas_evaluator.hh"
#include "vapor_pressure_water.hh"

#include "viscosity_evaluator.hh"
#include "viscosity_constant.hh"
#include "viscosity_water.hh"

#include "carbon_decomposition_rate_evaluator.hh"
#include "transport_decay_rate_evaluator.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<Evaluator, EOSEvaluator> EOSEvaluator::fac_("eos");

Utils::RegisteredFactory<EOS, EOSConstant> EOSConstant::factory_("constant");
Utils::RegisteredFactory<EOS, EOSIce> EOSIce::factory_("ice");
Utils::RegisteredFactory<EOS, EOSIdealGas> EOSIdealGas::factory_("ideal gas");
Utils::RegisteredFactory<EOS, EOSLinear> EOSLinear::factory_("linear");
Utils::RegisteredFactory<EOS, EOS_SW> EOS_SW::factory_("salt water");
Utils::RegisteredFactory<EOS, EOSVaporInGas> EOSVaporInGas::factory_("vapor in gas");
Utils::RegisteredFactory<EOS, EOSWater> EOSWater::factory_("liquid water");

Utils::RegisteredFactory<Evaluator, MolarFractionGasEvaluator> MolarFractionGasEvaluator::factory_(
  "molar fraction gas");
Utils::RegisteredFactory<VaporPressureRelation, VaporPressureWater> VaporPressureWater::factory_(
  "water vapor over water/ice");

Utils::RegisteredFactory<Evaluator, ViscosityEvaluator> ViscosityEvaluator::factory_("viscosity");
Utils::RegisteredFactory<ViscosityRelation, ViscosityConstant> ViscosityConstant::factory_(
  "constant");
Utils::RegisteredFactory<ViscosityRelation, ViscosityWater> ViscosityWater::factory_(
  "liquid water");

Utils::RegisteredFactory<Evaluator, CarbonDecomposeRateEvaluator>
  CarbonDecomposeRateEvaluator::reg_("carbon decomposition rate");
Utils::RegisteredFactory<Evaluator, TransportDecayRateEvaluator>
  TransportDecayRateEvaluator::reg_("transport decay rate");

} // namespace Relations
} // namespace Amanzi
