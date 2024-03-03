/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

//! Evaluates potential evapotranpiration (PET) using Priestley & Taylor formulation.
/*!

This implementation is based on models provided in the PRMS-IV, Version 4, see
pages 90-93, Equations 1-57 to 1-60

Requires the following dependencies:

.. _pet-priestley-taylor-evaluator-spec:
.. admonition:: pet-priestley-taylor-evaluator-spec:

   * `"include limiter`" ``[bool]`` **false** If true, multiply potential ET by
     a limiter to get an actual ET.
   * `"limiter number of dofs`" ``[int]`` **1** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides how many.
   * `"limiter dof`" ``[int]`` **0** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides which one to use.
   * `"include 1 - limiter`" ``[bool]`` **false** If true, multiply potential
     ET by 1 - a limiter (e.g. a limiter that partitions between two pools) to
     get actual ET.
   * `"1 - limiter number of dofs`" ``[int]`` **1** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides how many.
   * `"1 - limiter dof`" ``[int]`` **0** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides which one to use.
   * `"sublimate snow`" ``[bool]`` **false** If true, use latent heat of
      vaporization of snow, not water.

   KEYS:

   - `"air temperature`" **DOMAIN-air_temperature** Air temp, in [K]
   - `"surface temperature`" **DOMAIN-temperature** Ground or leaf temp, in [K].  Note this may be the
      same as air temperature.
   - `"elevation`" **DOMAIN-elevation** Elevation [m]
   - `"net radiation`" **DOMAIN-net_radiation** [W m^-2] Net radiation balance, positive to the ground.
   - `"limiter`" [-] See `"include limiter`" above.
   - `"1 - limiter`" [-] See `"include 1 - limiter`" above.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

namespace PriestleyTaylor {

KOKKOS_INLINE_FUNCTION
double
saturatedVaporPressure(double temp)
{
  // Westermann (2016), August–Roche–Magnus formula, Pa
  double tempC = temp - 273.15;
  if (tempC > 0) {
    return 611 * std::exp(17.62 * tempC / (tempC + 243.12));
  } else {
    return 611 * std::exp(22.46 * tempC / (tempC + 272.62));
  }
}


//
// PRMS-IV eqn 1-59, calculates ground heat flux density in units of
// [MJ m^-2 d^-1] given a daily-averaged ground and air temperature
// (in C or K).  We convert to W/m^2
//
KOKKOS_INLINE_FUNCTION
double
groundHeatFlux(double temp_ground, double temp_air)
{
  double G = -4.2 * (temp_ground - temp_air); // MJ/m^2/d
  return G * 1e6 / 86400;                     // W/m^2
}

//
// PRMS-IV eqn 1-58, calculates the slope of vapor pressure as a function of
// daily averaged air temperature [K], in [KPa C^-1]
//
KOKKOS_INLINE_FUNCTION
double
vaporPressureSlope(double temp_air)
{
  double tempC = temp_air - 273.15;                            // C
  double vp_sat = saturatedVaporPressure(temp_air); // Pa
  return 4098 * vp_sat / std::pow(tempC + 237.3, 2);           // Pa/C
}


//
// PRMS-IV eqn 1-57, calculates the psychrometric constant in [KPa C^-1] as a
// function of an elevation (lapse rate fixed) and a latent heat of
// vaporization in [cal gm^-1].
//
KOKKOS_INLINE_FUNCTION
double
psychrometricConstant(double lh_vap, double elev)
{
  double elev_ft = elev * 3.281;              // ft
  double lh_vap_calg = lh_vap / 1000 / 4.184; // cal/g
  double psy_kpaC = 1.6286 * (101.3 - (0.003215 * elev_ft)) / lh_vap_calg;
  // Kpa/C
  return psy_kpaC * 1000; // Pa/C
}


//
// PRMS-IV eqn 1-51, calculates the latent heat of vaporization [cal g^-1] as a
// function of the daily averaged air temperature [K] for liquid water
//
KOKKOS_INLINE_FUNCTION
double
latentHeatVaporization_water(double temp_air)
{
  double temp_f = 1.8 * (temp_air - 273.15) + 32; // F
  double lh_vap_calg = 597.3 - (0.5653 * temp_f); // cal/g
  return lh_vap_calg * 4.184 * 1000;              // J/kg
}


//
// PRMS-IV eqn 1-51, calculates the latent heat of vaporization [cal g^-1] as a
// function of the daily averaged air temperature [K] for snow -- note this is
// currently the same as the water value, but should get modified for snow!
//
double
latentHeatVaporization_snow(double temp_air)
{
  return latentHeatVaporization_water(temp_air); // J/kg
}



} // namespace PriestleyTaylor


class PETPriestleyTaylorEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit PETPriestleyTaylorEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  PETPriestleyTaylorEvaluator(const PETPriestleyTaylorEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new PETPriestleyTaylorEvaluator(*this));
  }

  std::string getType() const override { return eval_type; }

protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  static const std::string eval_type;

 protected:
  Key domain_;
  Key evap_type_;
  Key air_temp_key_;
  Key surf_temp_key_;
  Key elev_key_;
  Key rad_key_;
  Key limiter_key_;
  Key one_minus_limiter_key_;

  double pt_alpha_;
  bool limiter_, one_minus_limiter_;
  int limiter_nvecs_, one_minus_limiter_nvecs_;
  int limiter_dof_, one_minus_limiter_dof_;
  bool compatible_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, PETPriestleyTaylorEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
