/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//  Functions for calculating the snow / surface energy balance.


#include <iostream>
#include <cmath>
#include <algorithm>

#include "errors.hh"
#include "dbc.hh"
#include "Brent.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {

#define SWE_EPS 1.e-12
#define ENERGY_BALANCE_TOL 1.e-8


double
CalcAlbedoSnow(double density_snow)
{
  double AlSnow;
  if (density_snow <= 432.23309912785146) {
    AlSnow = 1.0 - 0.247 * std::pow(0.16 + 110 * std::pow(density_snow / 1000, 4), 0.5);
  } else {
    AlSnow = 0.6 - density_snow / 4600;
  }
  return AlSnow;
}

double
CalcRoughnessFactor(double snow_height, double Z_rough_bare, double Z_rough_snow)
{
  double Zfraction = snow_height <= 0.           ? 0. :
                     snow_height >= Z_rough_bare ? 1. :
                                                   1 - (Z_rough_bare - snow_height) / Z_rough_bare;
  return Z_rough_snow * Zfraction + Z_rough_bare * (1 - Zfraction);
}


std::pair<double, double>
IncomingRadiation(const MetData& met, double albedo)
{
  // Calculate incoming short-wave radiation
  double fQswIn = (1 - albedo) * met.QswIn;
  // Calculate incoming long-wave radiation
  double fQlwIn = met.QlwIn;
  return std::make_pair(fQswIn, fQlwIn);
}


//
// Calculate longwave from air temp and vapor pressure
// ------------------------------------------------------------------------------------------
double
IncomingLongwaveRadiation(double air_temp, double vapor_pressure_air)
{
  double vp_air_hPa = vapor_pressure_air / 100;
  double e_air = std::pow(vp_air_hPa, air_temp / 2016.);
  e_air = 1.08 * (1 - std::exp(-e_air));
  double longwave = e_air * c_stephan_boltzmann * std::pow(air_temp, 4);
  AMANZI_ASSERT(longwave > 0.);
  return longwave;
}


double
OutgoingLongwaveRadiation(double temp, double emissivity)
{
  // Calculate outgoing long-wave radiation
  return emissivity * c_stephan_boltzmann * std::pow(temp, 4);
}

double
BeersLawAbsorptivity(double k_extinction, double lai)
{
  return 1 - std::exp(-k_extinction * lai);
}

double
WindFactor(double Us, double Z_Us, double Z_rough, double KB)
{
  // Calculate D_h, D_e,
  return std::pow(c_von_Karman, 2) * Us /
         (std::pow(std::log(Z_Us / Z_rough), 2) + std::log(Z_Us / Z_rough) * KB);
}


double
StabilityFunction(double air_temp, double skin_temp, double Us, double Z_Us, double c_gravity)
{
  double Ri = c_gravity * Z_Us * (air_temp - skin_temp) / (air_temp * std::pow(Us, 2));
  if (Ri >= 0.) {
    // stable condition
    return 1. / (1 + 10 * Ri);
  } else {
    // Unstable condition
    return (1 - 10 * Ri);
  }
}


double
SaturatedVaporPressure(double temp)
{
  // Westermann (2016), August–Roche–Magnus formula, Pa
  double tempC = temp - 273.15;
  if (tempC > 0) {
    return 611 * std::exp(17.62 * tempC / (tempC + 243.12));
  } else {
    return 611 * std::exp(22.46 * tempC / (tempC + 272.62));
  }
}

double
SaturatedVaporPressureELM(double temp)
{
  // Saturated vapor pressure in [kPa] from CLM technical note
  double T = temp - 273.15;
  double coef_w[9] = { 6.1123516,     5.03109514e-1,  1.88369801e-2,  4.20547422e-4, 6.14396778e-6,
                       6.02780717e-8, 3.87940929e-10, 1.49436277e-12, 2.62655803e-15 };
  double coef_i[9] = { 6.11213467,    4.44007856e-1,  1.43064234e-2,   2.64461437e-4, 3.05903558e-6,
                       1.96237241e-8, 8.92344772e-11, -3.73208410e-13, 2.09339997e-16 };
  double* coef;
  if (T >= 0) coef = coef_w;
  else coef = coef_i;

  double res = coef[0];
  double Tn = T;
  for (int i = 1; i != 9; ++i) {
    res += coef[i] * Tn;
    Tn *= T;
  }
  // convert to Pa
  return 1e3 * res;
}

double
SaturatedSpecificHumidityELM(double temp, const ModelParams& params)
{
  double vp_sat = SaturatedVaporPressureELM(temp);
  return 0.622 * vp_sat / (params.P_atm - 0.378 * vp_sat);
}


double
VaporPressureGround(const GroundProperties& surf, const ModelParams& params)
{
  // Ho & Webb 2006
  double relative_humidity = -1;
  if (surf.pressure < params.P_atm) {
    // vapor pressure lowering
    double pc = params.P_atm - surf.pressure;
    relative_humidity = std::exp(-pc / (surf.density_w * c_R_ideal_gas * surf.temp));
  } else {
    relative_humidity = 1.;
  }
  return relative_humidity * SaturatedVaporPressure(surf.temp);
}


double
EvaporativeResistanceGround(const GroundProperties& surf,
                            const MetData& met,
                            double vapor_pressure_ground)
{
  // calculate evaporation prefactors
  if (met.vp_air > vapor_pressure_ground) { // condensation
    return 0.;
  } else {
    return surf.rsoil;
  }
}


double
SensibleHeat(double resistance_coef,
             double density_air,
             double Cp_air,
             double air_temp,
             double skin_temp)
{
  return resistance_coef * density_air * Cp_air * (air_temp - skin_temp);
}


double
LatentHeat(double resistance_coef,
           double density_air,
           double latent_heat_fusion,
           double vapor_pressure_air,
           double vapor_pressure_skin,
           double p_atm)
{
  return resistance_coef * density_air * latent_heat_fusion * 0.622 *
         (vapor_pressure_air - vapor_pressure_skin) / p_atm;
}

double
ConductedHeatIfSnow(double ground_temp, const SnowProperties& snow, const ModelParams& params)
{
  // Calculate heat conducted to ground, if snow
  double density = snow.density;
  if (density > 150) {
    // adjust for frost hoar
    density = 1. / ((0.90 / density) + (0.10 / 150));
  }
  double Ks = params.thermalK_freshsnow *
              std::pow(density / params.density_freshsnow, params.thermalK_snow_exp);
  return Ks * (snow.temp - ground_temp) / snow.height;
}


void
UpdateEnergyBalanceWithSnow_Inner(const GroundProperties& surf,
                                  const SnowProperties& snow,
                                  const MetData& met,
                                  const ModelParams& params,
                                  EnergyBalance& eb)
{
  // incoming radiation -- DONE IN OUTER

  // outgoing radiation
  eb.fQlwOut = OutgoingLongwaveRadiation(snow.temp, snow.emissivity);

  // sensible heat
  double Dhe = WindFactor(
    met.Us, met.Z_Us, CalcRoughnessFactor(snow.height, surf.roughness, snow.roughness), params.KB);
  double Sqig = StabilityFunction(met.air_temp, snow.temp, met.Us, met.Z_Us, params.gravity);
  eb.fQh = SensibleHeat(Dhe * Sqig, params.density_air, params.Cp_air, met.air_temp, snow.temp);

  // latent heat
  double vapor_pressure_skin = SaturatedVaporPressure(snow.temp);
  // KB here is to formulize the log ration between momentum roughness length and vapor roughness length
  // to add the effect of surface microtopography. Details see Gao et al. 2021 (WRR) Eq.(17) and similar
  // studies by Mölder & Lindroth 2001 (Agricultural and Forest Meteorology) but in heat transport.
  // Theories see Chapter 4 of book [Brutsaert, W. (1982). Evaporation into the atmosphere: Theory,
  // history, and applications]. Da0_a, Da0_b, Cd0_c, Cd0_d here are four fitting coefficients in the
  // formulization. There are several different values proposed by different studies (see Brutsaert 1982).
  // In Gao et al. (2021), these four coefficients are determined by lab experiments.
  // If set Da0_a, Cd0_c, Cd0_d to 0, it means that vapor roughness length is assumed equal to momentum
  // roughness length. Currently, this approach still needs testing or other formulization approaches
  // may also be added later for testing. So keep Da0_a = Cd0_c = Cd0_d = 0 for users.
  double u_star =
    met.Us * c_von_Karman /
    std::log(met.Z_Us / CalcRoughnessFactor(snow.height, surf.roughness, snow.roughness));
  double Re0 = params.density_air * u_star *
               CalcRoughnessFactor(snow.height, surf.roughness, snow.roughness) /
               params.dynamic_viscosity_air;
  double KB =
    (params.Da0_a * std::pow(Re0, params.Da0_b) - (params.Cd0_c * std::log(Re0) + params.Cd0_d)) *
    c_von_Karman;
  double Dhe_latent = WindFactor(
    met.Us, met.Z_Us, CalcRoughnessFactor(snow.height, surf.roughness, snow.roughness), KB);
  double coef = std::min(Dhe_latent * Sqig, 1.0);
  eb.fQe = LatentHeat(
    coef, params.density_air, params.H_sublimation, met.vp_air, vapor_pressure_skin, params.P_atm);


  // conducted heat
  eb.fQc = ConductedHeatIfSnow(surf.temp, snow, params);

  // balance of energy goes into melting
  eb.fQm = eb.fQswIn + eb.fQlwIn - eb.fQlwOut + eb.fQh - eb.fQc + eb.fQe;
}

EnergyBalance
UpdateEnergyBalanceWithSnow(const GroundProperties& surf,
                            const MetData& met,
                            const ModelParams& params,
                            SnowProperties& snow)
{
  EnergyBalance eb;

  // snow on the ground, solve for snow temperature
  std::tie(eb.fQswIn, eb.fQlwIn) = IncomingRadiation(met, snow.albedo);
  snow.temp = DetermineSnowTemperature(surf, met, params, snow, eb);

  if (snow.temp > 273.15) {
    // limit snow temp to 0, then melt with the remaining energy
    snow.temp = 273.15;
    UpdateEnergyBalanceWithSnow_Inner(surf, snow, met, params, eb);
    eb.error = 0.;
  } else {
    // snow not melting
    UpdateEnergyBalanceWithSnow_Inner(surf, snow, met, params, eb);
    eb.error = eb.fQm;
    eb.fQm = 0.;
  }

  return eb;
}


EnergyBalance
UpdateEnergyBalanceWithoutSnow(const GroundProperties& surf,
                               const MetData& met,
                               const ModelParams& params)
{
  EnergyBalance eb;

  // incoming radiation
  std::tie(eb.fQswIn, eb.fQlwIn) = IncomingRadiation(met, surf.albedo);

  // outgoing radiation
  eb.fQlwOut = OutgoingLongwaveRadiation(surf.temp, surf.emissivity);

  // potentially have incoming precip to melt
  if (surf.temp > 273.65) {
    eb.fQm = (met.Ps + surf.snow_death_rate) * surf.density_w * params.H_fusion;
  } else if (surf.temp <= 273.15) {
    eb.fQm = 0.;
  } else {
    double Em = (met.Ps + surf.snow_death_rate) * surf.density_w * params.H_fusion;
    eb.fQm = Em * (surf.temp - 273.15) / (0.5);
  }

  // sensible heat
  double Dhe = WindFactor(met.Us, met.Z_Us, surf.roughness, params.KB);
  double Sqig = StabilityFunction(met.air_temp, surf.temp, met.Us, met.Z_Us, params.gravity);
  eb.fQh = SensibleHeat(Dhe * Sqig, params.density_air, params.Cp_air, met.air_temp, surf.temp);

  // latent heat
  double vapor_pressure_skin = VaporPressureGround(surf, params);
  double u_star = met.Us * c_von_Karman / std::log(met.Z_Us / surf.roughness);
  double Re0 = params.density_air * u_star * surf.roughness / params.dynamic_viscosity_air;
  double KB =
    (params.Da0_a * std::pow(Re0, params.Da0_b) - (params.Cd0_c * std::log(Re0) + params.Cd0_d)) *
    c_von_Karman;
  double Dhe_latent = WindFactor(met.Us, met.Z_Us, surf.roughness, KB);
  double Rsoil = EvaporativeResistanceGround(surf, met, vapor_pressure_skin);
  double coef = std::min(1.0 / (Rsoil + 1.0 / (Dhe_latent * Sqig)), 1.0);

  // positive is condensation
  eb.fQe = LatentHeat(coef,
                      params.density_air,
                      surf.unfrozen_fraction * params.H_vaporization +
                        (1 - surf.unfrozen_fraction) * params.H_sublimation,
                      met.vp_air,
                      vapor_pressure_skin,
                      params.P_atm);

  // fQc is the energy balance -- this is not used in this branch (e.g. no snow
  // branch) but is kept here anyway for use by diagnostic variables
  eb.fQc = eb.fQswIn + eb.fQlwIn - eb.fQlwOut + eb.fQh + eb.fQe;
  return eb;
}

// Snow temperature calculation.
double
DetermineSnowTemperature(const GroundProperties& surf,
                         const MetData& met,
                         const ModelParams& params,
                         SnowProperties& snow,
                         EnergyBalance& eb,
                         std::string method)
{
  SnowTemperatureFunctor_ func(&surf, &snow, &met, &params, &eb);
  double left, right;
  double res_left, res_right;

  double res_init = func(surf.temp);
  if (res_init < 0.) {
    right = surf.temp;
    res_right = res_init;

    left = surf.temp - 1.;
    res_left = func(left);
    while (res_left < 0.) {
      right = left;
      res_right = res_left;
      left = left - 1.;
      res_left = func(left);
    }
  } else {
    left = surf.temp;
    res_left = res_init;

    right = surf.temp + 1.;
    res_right = func(right);
    while (res_right > 0.) {
      left = right;
      res_left = res_right;
      right = right + 1.;
      res_right = func(right);
    }
  }

  int my_max_it = 100;
  int max_it(my_max_it);
  double result(0.);
  if (method == "bisection") {
    Errors::Message msg("SurfaceEnergyBalance: root finding method \"bisection\" is not longer "
                        "supported -- use \"brent\"");
    Exceptions::amanzi_throw(msg);
  } else if (method == "toms") {
    Errors::Message msg("SurfaceEnergyBalance: root finding method \"bisection\" is not longer "
                        "supported -- use \"brent\"");
    Exceptions::amanzi_throw(msg);
  } else if (method == "brent") {
    result = Utils::findRootBrent(func, left, right, ENERGY_BALANCE_TOL, &max_it);
  } else {
    Errors::Message emsg;
    emsg << "SurfaceEnergyBalance: invalid solver method \"" << method << "\", use \"brent\"";
    Exceptions::amanzi_throw(emsg);
  }

  if (max_it >= my_max_it) throw("Nonconverged Surface Energy Balance");
  // Call the function again to set the fluxes.
  func(result);
  return result;
}


MassBalance
UpdateMassBalanceWithSnow(const GroundProperties& surf,
                          const ModelParams& params,
                          const EnergyBalance& eb)
{
  MassBalance mb;

  // Melt rate given by available energy rate divided by heat of fusion.
  mb.Mm = eb.fQm / (surf.density_w * params.H_fusion);

  // Snow balance
  mb.Me = eb.fQe / (surf.density_w * params.H_sublimation);
  return mb;
}

MassBalance
UpdateMassBalanceWithoutSnow(const GroundProperties& surf,
                             const ModelParams& params,
                             const EnergyBalance& eb)
{
  MassBalance mb;
  mb.Mm = eb.fQm / (surf.density_w * params.H_fusion);
  mb.Me = eb.fQe / (surf.density_w * (surf.unfrozen_fraction * params.H_vaporization +
                                      (1 - surf.unfrozen_fraction) * params.H_sublimation));
  return mb;
}

FluxBalance
UpdateFluxesWithoutSnow(const GroundProperties& surf,
                        const MetData& met,
                        const ModelParams& params,
                        const EnergyBalance& eb,
                        const MassBalance& mb,
                        bool model_1p1)
{
  FluxBalance flux;

  // mass to surface is precip and melting first
  // partition all fluxes?
  // double water_flux = met.Pr + mb.Mm + mb.Me;
  // flux.M_surf = 0.;

  // or just partition evaporation?
  double water_flux = mb.Me;
  flux.M_surf = met.Pr + mb.Mm;

  // Energy to surface.
  double Train = std::max(0., met.air_temp - 273.15);
  flux.E_surf = eb.fQswIn + eb.fQlwIn - eb.fQlwOut + eb.fQh        // purely energy fluxes
              - eb.fQm                                             // energy put into melting snow
              + surf.density_w * met.Pr * Train * params.Cv_water; // energy advected in by rainfall

  // zero subsurf values
  flux.M_subsurf = 0.;
  flux.E_subsurf = 0.;

  // allocate water_flux to surface or subsurface
  double evap_to_subsurface_fraction;
  if (model_1p1) {
    // NOTE: this old model allows indepedent values of evap_transition_width
    // (which governs where the flux goes, in units of Pa) and
    // water_transition_depth (which governs where the
    // saturation/pressure are used to calculate how much water to take, in
    // units of [m]).  This was deprecated because it allowed inconsistencies.
    if (mb.Me < 0) {
      if (surf.pressure >= params.P_atm + params.evap_transition_width) {
        evap_to_subsurface_fraction = 0.;
      } else if (surf.pressure < params.P_atm) {
        evap_to_subsurface_fraction = 1.;
      } else {
        evap_to_subsurface_fraction =
          (params.P_atm + params.evap_transition_width - surf.pressure) /
          (params.evap_transition_width);
      }
    } else {
      evap_to_subsurface_fraction = 0.;
    }
  } else {
    evap_to_subsurface_fraction =
      1 - std::min(1.0, surf.ponded_depth / surf.water_transition_depth);
  }
  AMANZI_ASSERT(evap_to_subsurface_fraction >= 0. && evap_to_subsurface_fraction <= 1.);

  flux.M_surf += (1 - evap_to_subsurface_fraction) * water_flux;
  flux.M_subsurf += evap_to_subsurface_fraction * water_flux;

  // enthalpy of evap/condensation always goes entirely to surface
  //
  // This is fine because diffusion of energy always works, and doing otherwise
  // sets up local minima for energy in the top cell, which breaks code.
  flux.E_surf += eb.fQe; // v1.0

  // // enthalpy of evap/condensation -- v0.88
  // flux.E_surf += (1-evap_to_subsurface_fraction) * eb.fQe;
  // flux.E_subsurf += evap_to_subsurface_fraction * eb.fQe;

  // snow mass change
  flux.M_snow = met.Ps - mb.Mm;
  return flux;
}


FluxBalance
UpdateFluxesWithSnow(const GroundProperties& surf,
                     const MetData& met,
                     const ModelParams& params,
                     const SnowProperties& snow,
                     const EnergyBalance& eb,
                     const MassBalance& mb)
{
  FluxBalance flux;

  // mass to surface is precip and evaporation
  flux.M_surf = met.Pr + mb.Mm;
  if (mb.Mm > 0.) AMANZI_ASSERT(snow.density > 99.);
  flux.M_snow = met.Ps + mb.Me - mb.Mm;

  // Energy to surface.
  double Train = std::max(0., met.air_temp - 273.15);
  flux.E_surf = eb.fQc                                             // conducted to ground
              + surf.density_w * met.Pr * Train * params.Cv_water; // rain enthalpy
  // + 0 // enthalpy of meltwater at 0C.
  return flux;
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
