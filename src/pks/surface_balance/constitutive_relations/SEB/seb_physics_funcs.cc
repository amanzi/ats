/*
  Functions for calculating the snow / surface energy balance.
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include "boost/math/tools/roots.hpp"

#include "dbc.hh"

#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {

#define SWE_EPS 1.e-12
#define ENERGY_BALANCE_TOL 1.e-8


double CalcAlbedoSnow(double density_snow) {
  double AlSnow;
  if (density_snow <= 432.23309912785146) {
    AlSnow = 1.0 - 0.247 * std::pow(0.16 + 110*std::pow(density_snow/1000, 4), 0.5);
  } else {
    AlSnow = 0.6 - density_snow / 4600;
  }
  return AlSnow;
}

double CalcRoughnessFactor(double snow_height, double Z_rough_bare, double Z_rough_snow)
{
  double Zfraction = snow_height <= 0. ? 0. : snow_height >= Z_rough_bare ? 1. :
                     1 - (Z_rough_bare-snow_height)/Z_rough_bare;
  return Z_rough_snow * Zfraction + Z_rough_bare * (1-Zfraction);
}


std::pair<double,double> IncomingRadiation(const MetData& met, double albedo)
{
  // Calculate incoming short-wave radiation
  double fQswIn = (1 - albedo) * met.QswIn;

  // Calculate incoming long-wave radiation
  double fQlwIn = met.QlwIn;

  return std::make_pair(fQswIn, fQlwIn);
}


//
// Calculate longwave from air temp and relative humidity
// ------------------------------------------------------------------------------------------
double CalcIncomingLongwave(double air_temp, double relative_humidity, double c_stephan_boltzmann) {
  double e_air = std::pow(10 * VaporPressureAir(air_temp, relative_humidity), air_temp / 2016.);
  e_air = 1.08 * (1 - std::exp(-e_air));
  double longwave = e_air * c_stephan_boltzmann * std::pow(air_temp,4);
  AMANZI_ASSERT(longwave > 0.);
  return longwave;
}


double OutgoingRadiation(double temp, double emissivity, double c_stephan_boltzmann)
{
  // Calculate outgoing long-wave radiation
  return emissivity * c_stephan_boltzmann * std::pow(temp,4);
}

double WindFactor(double Us, double Z_Us, double Z_rough, double c_von_Karman)
{
  // Calculate D_h, D_e, 
  return std::pow(c_von_Karman,2) * Us / std::pow(std::log(Z_Us / Z_rough), 2);
}


double StabilityFunction(double air_temp, double skin_temp, double Us,
                         double Z_Us, double c_gravity)
{
  double Ri  = c_gravity * Z_Us * (air_temp - skin_temp) / (air_temp * std::pow(Us,2));
  if (Ri >= 0.) {
    // stable condition
    return 1. / (1 + 10*Ri);
  } else {
    // Unstable condition
    return (1 - 10*Ri);
  }
}


double SaturatedVaporPressure(double temp)
{
  // Sat vap. press o/water Dingman D-7 (Bolton, 1980)
  // *** (Bolton, 1980) Calculates vapor pressure in [kPa]  ****
  double tempC = temp - 273.15;
  return 0.6112 * std::exp(17.67 * tempC / (tempC + 243.5));
}

double VaporPressureAir(double air_temp, double relative_humidity)
{
  return SaturatedVaporPressure(air_temp) * relative_humidity;
}

double VaporPressureGround(const GroundProperties& surf, const ModelParams& params)
{
  // Ho & Webb 2006
  double relative_humidity = -1;
  if (surf.pressure < params.Apa * 1000.) {
    // vapor pressure lowering
    double pc = 1000.*params.Apa - surf.pressure;
    relative_humidity = std::exp(-pc / (surf.density_w * params.R_ideal_gas * surf.temp));
  } else {
    relative_humidity = 1.;
  }
  return relative_humidity * SaturatedVaporPressure(surf.temp);
}


double EvaporativeResistanceGround(const GroundProperties& surf,
        const MetData& met,
        const ModelParams& params, 
        double vapor_pressure_air, double vapor_pressure_ground)
{
  // calculate evaporation prefactors
  if (vapor_pressure_air > vapor_pressure_ground) { // condensation
    return 0.;
  } else {
    return EvaporativeResistanceCoef(surf.saturation_gas,
            surf.porosity, surf.dz, params.Clapp_Horn_b);
  }
}

double EvaporativeResistanceCoef(double saturation_gas,
        double porosity, double dessicated_zone_thickness, double Clapp_Horn_b) {
  double Rsoil;
  if (saturation_gas == 0.) {
    Rsoil = 0.; // ponded water
  } else {
    // Equation for reduced vapor diffusivity
    // See Sakagucki and Zeng 2009 eqaution (9) and Moldrup et al., 2004. 
    double vp_diffusion = 0.000022 * (std::pow(porosity,2))
                          * std::pow((1-(0.0556/porosity)),(2+3*Clapp_Horn_b));
    // Sakagucki and Zeng 2009 eqaution (10)
    double L_Rsoil = std::exp(std::pow(saturation_gas, 5));
    L_Rsoil = dessicated_zone_thickness * (L_Rsoil -1) * (1/(std::exp(1.)-1));
    Rsoil = L_Rsoil/vp_diffusion;
  }
  return Rsoil;
}


double SensibleHeat(double resistance_coef,
                    double density_air,
                    double Cp_air,
                    double air_temp,
                    double skin_temp) {
  return resistance_coef * density_air * Cp_air * (air_temp - skin_temp);
}


double LatentHeat(double resistance_coef,
                  double density_air,  /// this should be w?
                  double latent_heat_fusion,
                  double vapor_pressure_air,
                  double vapor_pressure_skin,
                  double Apa) {
  AMANZI_ASSERT(resistance_coef <= 1.);
  return resistance_coef * density_air * latent_heat_fusion * 0.622
      * (vapor_pressure_air - vapor_pressure_skin) / Apa;
}

double ConductedHeatIfSnow(double ground_temp,
                           const SnowProperties& snow, const ModelParams& params)
{
  // Calculate heat conducted to ground, if snow
  double density = snow.density;
  if (density > 150) {
    // adjust for frost hoar
    density = 1. / ((0.90/density) + (0.10/150));
  }
  double Ks = params.thermalK_freshsnow * std::pow(density/params.density_freshsnow, params.thermalK_snow_exp);
  return Ks * (snow.temp - ground_temp) / snow.height;
}


void UpdateEnergyBalanceWithSnow_Inner(const GroundProperties& surf,
        const SnowProperties& snow,
        const MetData& met,
        const ModelParams& params,
        EnergyBalance& eb)
{
  // incoming radiation -- DONE IN OUTER

  // outgoing radiation
  eb.fQlwOut = OutgoingRadiation(snow.temp, snow.emissivity, params.stephB);

  // sensible heat
  double Dhe = WindFactor(met.Us, met.Z_Us, CalcRoughnessFactor(snow.height, surf.roughness, snow.roughness), params.VKc);
  double Sqig = StabilityFunction(met.air_temp, snow.temp, met.Us, met.Z_Us, params.gravity);
  eb.fQh = SensibleHeat(Dhe * Sqig, params.density_air, params.Cp_air, met.air_temp, snow.temp);

  // latent heat
  double vapor_pressure_air = VaporPressureAir(met.air_temp, met.relative_humidity);
  double vapor_pressure_skin = SaturatedVaporPressure(snow.temp);
  eb.fQe = LatentHeat(Dhe * Sqig, params.density_air, params.Ls, vapor_pressure_air, vapor_pressure_skin,
                      params.Apa);

  // conducted heat
  eb.fQc = ConductedHeatIfSnow(surf.temp, snow, params);

  // balance of energy goes into melting
  eb.fQm = eb.fQswIn + eb.fQlwIn - eb.fQlwOut + eb.fQh - eb.fQc + eb.fQe;
}

EnergyBalance UpdateEnergyBalanceWithSnow(const GroundProperties& surf,
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


EnergyBalance UpdateEnergyBalanceWithoutSnow(const GroundProperties& surf,
        const MetData& met,
        const ModelParams& params)
{
  EnergyBalance eb;

  // incoming radiation
  std::tie(eb.fQswIn, eb.fQlwIn) = IncomingRadiation(met, surf.albedo);

  // outgoing radiation
  eb.fQlwOut = OutgoingRadiation(surf.temp, surf.emissivity, params.stephB);

  // potentially have incoming precip to melt
  if (surf.temp > 273.65) {
    eb.fQm = (met.Ps + surf.snow_death_rate) * surf.density_w * params.Hf;
  } else if (surf.temp <= 273.15) {
    eb.fQm = 0.;
  } else {
    double Em = (met.Ps + surf.snow_death_rate) * surf.density_w * params.Hf;
    eb.fQm = Em * (surf.temp - 273.15) / (0.5);
  }
  
  // sensible heat
  double Dhe = WindFactor(met.Us, met.Z_Us, surf.roughness, params.VKc);
  double Sqig = StabilityFunction(met.air_temp, surf.temp, met.Us, met.Z_Us, params.gravity);
  eb.fQh = SensibleHeat(Dhe*Sqig, params.density_air, params.Cp_air, met.air_temp, surf.temp);

  // latent heat
  double vapor_pressure_air = VaporPressureAir(met.air_temp, met.relative_humidity);
  double vapor_pressure_skin = VaporPressureGround(surf, params);

  double Rsoil = EvaporativeResistanceGround(surf, met, params, vapor_pressure_air, vapor_pressure_skin);
  double coef = 1.0 / (Rsoil + 1.0/(Dhe*Sqig));


  // NUMERICAL DOWNREGULATION due to difficulty evaporating ice...
  //  coef *= surf.unfrozen_fraction;
  
  // positive is condensation
  eb.fQe = LatentHeat(coef, params.density_air, 
		      surf.unfrozen_fraction * params.Le + (1-surf.unfrozen_fraction) * params.Ls,
		      vapor_pressure_air, vapor_pressure_skin, params.Apa);

  //This is a numerical fix to turnoff condensation when ground is not unfrozen.
  if (!params.condensation_frozen) {
    if (surf.unfrozen_fraction != 1 && eb.fQe > 0) eb.fQe = 0;
  }
  if (!params.evaporation_frozen) {
    if (surf.unfrozen_fraction < 1) eb.fQe *= surf.unfrozen_fraction;
  }
  // fQc is the energy conducted between surface and snow layers, but there is no snow here
  eb.fQc = 0.;
  return eb;
}

// Snow temperature calculation.
double DetermineSnowTemperature(const GroundProperties& surf,
        const MetData& met,
        const ModelParams& params, 
        SnowProperties& snow,
        EnergyBalance& eb,
        std::string method)
{
  SnowTemperatureFunctor_ func(&surf, &snow, &met, &params, &eb);
  Tol_ tol(ENERGY_BALANCE_TOL);
  boost::uintmax_t max_it(100);
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

  //  std::cout << "Determining snow temp in interval (" << left << "," << right << ")" << std::endl;
  
  std::pair<double,double> result;
  auto my_max_it = max_it;
  if (method == "bisection") {
    result = boost::math::tools::bisect(func, left, right, tol, max_it);
  } else if (method == "toms") {
    result = boost::math::tools::toms748_solve(func, left, right, res_left, res_right, tol, max_it);
  }

  if (max_it >= my_max_it) throw("Nonconverged Surface Energy Balance");
  // boost algorithms calculate the root such that the root is contained within the interval
  // [result.first, result.second], and that width of that interval is less than TOL.
  // We choose to set the solution as the center of that interval.
  // Call the function again to set the fluxes.
  double solution = (result.first + result.second)/2.;
  //  std::cout << "  Got T=" << solution << " with SW = " << eb.fQswIn << ", LW = " << eb.fQlwIn << ", LWout = " << eb.fQlwOut << ", Sens = " << eb.fQh << ", Lat = " << eb.fQe << ", Cond = " << eb.fQc << ", melt = " << eb.fQm << std::endl;

  return solution;
}


MassBalance UpdateMassBalanceWithSnow(const GroundProperties& surf,
        const ModelParams& params, const EnergyBalance& eb)
{
  MassBalance mb;

  // Melt rate given by available energy rate divided by heat of fusion.
  mb.Mm = eb.fQm / (surf.density_w * params.Hf);

  // Snow balance
  mb.Me = eb.fQe / (surf.density_w * params.Ls);
  return mb;
}

MassBalance UpdateMassBalanceWithoutSnow(const GroundProperties& surf,
        const ModelParams& params, const EnergyBalance& eb)
{
  MassBalance mb;
  mb.Mm = eb.fQm / (surf.density_w * params.Hf);
  mb.Me = eb.fQe / (surf.density_w * (surf.unfrozen_fraction * params.Le + (1-surf.unfrozen_fraction) * params.Ls));
  return mb;
}

FluxBalance UpdateFluxesWithoutSnow(const GroundProperties& surf,
        const MetData& met, const ModelParams& params, const EnergyBalance& eb, const MassBalance& mb)
{
  FluxBalance flux;

  // mass to surface is precip and melting first
  flux.M_surf = met.Pr + mb.Mm;

  // Energy to surface.
  double Train = std::max(0., met.air_temp - 273.15);
  flux.E_surf = eb.fQswIn + eb.fQlwIn - eb.fQlwOut + eb.fQh // purely energy fluxes
                - eb.fQm   // energy put into melting snow
                + surf.density_w * met.Pr * Train * params.Cv_water; // energy advected in by rainfall

  // zero subsurf values
  flux.M_subsurf = 0.;
  flux.E_subsurf = 0.;

  // At this point we have Mass and Energy fluxes but not including
  // evaporation, which we have to allocate to surface or subsurface.
  double evap_to_subsurface_fraction = 0.;
  if (mb.Me < 0) {
    if (surf.pressure >= 1000.*params.Apa + params.evap_transition_width) {
      evap_to_subsurface_fraction = 0.;
    } else if (surf.pressure < 1000.*params.Apa) {
      evap_to_subsurface_fraction = 1.;
    } else {
      evap_to_subsurface_fraction = (1000.*params.Apa + params.evap_transition_width - surf.pressure) / (params.evap_transition_width);
    }
  }
  AMANZI_ASSERT(evap_to_subsurface_fraction >= 0. && evap_to_subsurface_fraction <= 1.);
  flux.M_surf += (1. - evap_to_subsurface_fraction) * mb.Me;
  flux.M_subsurf += evap_to_subsurface_fraction * mb.Me;

  // enthalpy of evap/condensation always goes entirely to surface
  //
  // This is fine because diffusion of energy always works, and doing otherwise
  // sets up local minima for energy in the top cell, which breaks code.
  //flux.E_surf += eb.fQe; // v1.0

  // enthalpy of evap/condensation -- v0.88
  flux.E_surf += (1-evap_to_subsurface_fraction) * eb.fQe;
  flux.E_subsurf += evap_to_subsurface_fraction * eb.fQe;
  
  // snow mass change
  flux.M_snow = met.Ps - mb.Mm;
  return flux;
}


FluxBalance UpdateFluxesWithSnow(const GroundProperties& surf,
        const MetData& met, const ModelParams& params, const SnowProperties& snow,
        const EnergyBalance& eb, const MassBalance& mb)
{
  FluxBalance flux;

  // mass to surface is precip and evaporation
  flux.M_surf = met.Pr + mb.Mm;
  if (mb.Mm > 0.) AMANZI_ASSERT(snow.density > 99.);
  flux.M_snow = met.Ps + mb.Me - mb.Mm;

  // Energy to surface.
  double Train = std::max(0., met.air_temp - 273.15);
  flux.E_surf = eb.fQc   // conducted to ground  
                + surf.density_w * met.Pr * Train * params.Cv_water; // rain enthalpy
               // + 0 // enthalpy of meltwater at 0C.
  return flux;
}





} // namespace
} // namespace
} // namespace
