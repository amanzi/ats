/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Functions for calculating the snow / surface energy balance.
*/

#pragma once

#include <cmath>

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
  double tempC = temp_air - 273.15;                  // C
  double vp_sat = saturatedVaporPressure(temp_air);  // Pa
  return 4098 * vp_sat / std::pow(tempC + 237.3, 2); // Pa/C
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
KOKKOS_INLINE_FUNCTION
double
latentHeatVaporization_snow(double temp_air)
{
  return latentHeatVaporization_water(temp_air); // J/kg
}


} // namespace PriestleyTaylor

namespace Functions {

static const double c_stephan_boltzmann = 0.00000005670373; // Stephan-Boltzmann constant
static const double c_von_Karman = 0.41;                    // [-] Von Karman constant
static const double c_R_ideal_gas = 461.52;                 // [Pa m^3 kg^-1 K^-1]


KOKKOS_INLINE_FUNCTION
double
beersLawAbsorptivity(double k_extinction, double lai)
{
  return 1 - exp(-k_extinction * lai);
}

KOKKOS_INLINE_FUNCTION
double
outgoingLongwaveRadiation(double temp, double emissivity)
{
  // Calculate outgoing long-wave radiation using the Boltzmann equation
  return emissivity * c_stephan_boltzmann * std::pow(temp, 4);
}


KOKKOS_INLINE_FUNCTION
double
atmosphereLongwaveRadiation(double air_temp, double vapor_pressure_air)
{
  // computes the atmosphere's outgoing longwave
  double vp_air_hPa = vapor_pressure_air / 100;
  double e_air = std::pow(vp_air_hPa, air_temp / 2016.);
  e_air = 1.08 * (1 - exp(-e_air));
  return outgoingLongwaveRadiation(air_temp, e_air);
}


KOKKOS_INLINE_FUNCTION
double
albedoSnow(const double density_snow)
{
  double a_snow;
  if (density_snow <= 432.23309912785146) {
    a_snow = 1.0 - 0.247 * pow(0.16 + 110 * pow(density_snow / 1000, 4), 0.5);
  } else {
    a_snow = 0.6 - density_snow / 4600;
  }
  return a_snow;
}

} // namespace Functions
} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
