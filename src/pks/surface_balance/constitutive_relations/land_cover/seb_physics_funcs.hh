/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Functions for calculating the snow / surface energy balance.

This formulation is based on Ling & Zhang 2004 which is tuned slightly toward
Arctic/permafrost applications:

- Ling, F., & Zhang, T. (2004). A numerical model for surface energy balance
and thermal regime of the active layer and permafrost containing unfrozen
water. Cold Regions Science and Technology, 38(1), 1-15.

and is documented in Appendix B of Atchley et al 2015.

- Atchley, A. L., Painter, S. L., Harp, D. R., Coon, E. T., Wilson, C. J.,
  Liljedahl, A. K., & Romanovsky, V. E. (2015). Using field observations to
  inform thermal hydrology models of permafrost dynamics with ATS
  (v0. 83). Geoscientific Model Development, 8(9), 2701-2722.

For evaporative fluxes in particular, the model is relatively sensitive to
relative humidity and, to a lesser extent, wind speed (in forcing datasets) and
to the roughness lengths (in parameters).

*/

#pragma once

#include <cmath>
#include <string>

#include "VerboseObject.hh"
#include "seb_physics_defs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Main SEB functions

//
// Determine the albedo of snow as a function of density.
// ------------------------------------------------------------------------------------------
double
CalcAlbedoSnow(double density_snow);

//
// Determine the surface roughness
// ------------------------------------------------------------------------------------------
double
CalcRoughnessFactor(double snow_height, double Z_rough_bare, double Z_rough_snow);


//
// Calculate longwave from air temp and vapor pressure
// ------------------------------------------------------------------------------------------
double
IncomingLongwaveRadiation(double air_temp, double vapor_pressure_air);

//
// Calculates incoming shortwave and longwave radiation incident on surface
// ------------------------------------------------------------------------------------------
std::pair<double, double>
IncomingRadiation(const MetData& met, double albedo);

//
// Calculates outgoing longwave radiation
// ------------------------------------------------------------------------------------------
double
OutgoingLongwaveRadiation(double temp, double emissivity);

//
// Beer's law for radiation attenuation through a single-layer canopy
// ------------------------------------------------------------------------------------------
double
BeersLawAbsorptivity(double k_extinction, double lai);

//
// Wind speed term D_he
// ------------------------------------------------------------------------------------------
double
WindFactor(double Us, double Z_Us, double Z_rough, double KB);

//
// Stability of convective overturning term Zeta AKA Sqig
// ------------------------------------------------------------------------------------------
double
StabilityFunction(double air_temp, double skin_temp, double Us, double Z_Us, double c_gravity);


//
// Westermann 2016, saturated vapor pressure over water/ice
// In [Pa]
// ------------------------------------------------------------------------------------------
double
SaturatedVaporPressure(double temp);
double
SaturatedVaporPressureELM(double temp);
double
SaturatedSpecificHumidityELM(double temp);


//
// Partial pressure of water vapor in gaseous phase, in the soil.
// After Ho & Webb 2006
// In [Pa]
// ------------------------------------------------------------------------------------------
double
VaporPressureGround(const GroundProperties& surf, const ModelParams& params);


//
// Diffusion of vapor pressure limiter on evaporation.
// After Sakagucki and Zeng 2009 eqaution (10)
// ------------------------------------------------------------------------------------------
double
EvaporativeResistanceGround(const GroundProperties& surf,
                            const MetData& met,
                            double vapor_pressure_ground);


//
// Basic sensible heat.
// ------------------------------------------------------------------------------------------
double
SensibleHeat(double resistance_coef,
             double density_air,
             double Cp_air,
             double air_temp,
             double skin_temp);

//
// Basic latent heat.
// ------------------------------------------------------------------------------------------
double
LatentHeat(double resistance_coef,
           double density_air, /// this should be w?
           double latent_heat_fusion,
           double vapor_pressure_air,
           double vapor_pressure_skin,
           double Apa);

//
// Heat conducted to ground via simple diffusion model between snow and skin surface.
// ------------------------------------------------------------------------------------------
double
ConductedHeatIfSnow(double ground_temp, const SnowProperties& snow);

//
// Update the energy balance, solving for the amount of heat available to melt snow.
//
// NOTE, this should not be used directly -- instead it is called within the loop solving for
// snow temperature.
// ------------------------------------------------------------------------------------------
void
UpdateEnergyBalanceWithSnow_Inner(const GroundProperties& surf,
                                  const SnowProperties& snow,
                                  const MetData& met,
                                  const ModelParams& params,
                                  EnergyBalance& eb);

//
// Determine the snow temperature by solving for energy balance, i.e. the snow
// temp at equilibrium.  Assumes no melting (and therefore T_snow calculated
// can be greater than 0 C.
// ------------------------------------------------------------------------------------------
double
DetermineSnowTemperature(const GroundProperties& surf,
                         const MetData& met,
                         const ModelParams& params,
                         SnowProperties& snow,
                         EnergyBalance& eb,
                         std::string method = "toms");


//
// Update the energy balance, solving for the amount of heat conducted to the ground.
//
// NOTE, this CAN be used directly.
// ------------------------------------------------------------------------------------------
EnergyBalance
UpdateEnergyBalanceWithSnow(const GroundProperties& surf,
                            const MetData& met,
                            const ModelParams& params,
                            SnowProperties& snow);

//
// Update the energy balance, solving for the amount of heat conducted to the ground.
//
// NOTE, this CAN be used directly.
// ------------------------------------------------------------------------------------------
EnergyBalance
UpdateEnergyBalanceWithoutSnow(const GroundProperties& surf,
                               const MetData& met,
                               const ModelParams& params);

//
// Given an energy balance, determine the resulting mass changes between
// precip, evaporation, melt, etc, with snow.
// ------------------------------------------------------------------------------------------
MassBalance
UpdateMassBalanceWithSnow(const GroundProperties& surf,
                          const ModelParams& params,
                          const EnergyBalance& eb);

//
// Given an energy balance, determine the resulting mass changes between
// precip, evaporation, melt, etc, with snow.
// ------------------------------------------------------------------------------------------
MassBalance
UpdateMassBalanceWithoutSnow(const GroundProperties& surf,
                             const ModelParams& params,
                             const EnergyBalance& eb);


//
// Given an energy balance and a mass balance, accumulate these into sources
// for surf and subsurf.
// ------------------------------------------------------------------------------------------
FluxBalance
UpdateFluxesWithSnow(const GroundProperties& surf,
                     const MetData& met,
                     const ModelParams& params,
                     const SnowProperties& snow,
                     const EnergyBalance& eb,
                     const MassBalance& mb);

//
// Given an energy balance and a mass balance, accumulate these into sources
// for surf and subsurf.
// ------------------------------------------------------------------------------------------
FluxBalance
UpdateFluxesWithoutSnow(const GroundProperties& surf,
                        const MetData& met,
                        const ModelParams& params,
                        const EnergyBalance& eb,
                        const MassBalance& mb,
                        bool model_1p1 = false);


// Calculation of a snow temperature requires a root-finding operation, for
// which we use a functor.
class SnowTemperatureFunctor_ {
 public:
  explicit SnowTemperatureFunctor_(GroundProperties const* const surf,
                                   SnowProperties* const snow,
                                   MetData const* const met,
                                   ModelParams const* const params,
                                   EnergyBalance* const eb)
    : surf_(surf), params_(params), snow_(snow), met_(met), eb_(eb)
  {}

  double operator()(double temp)
  {
    snow_->temp = temp;
    UpdateEnergyBalanceWithSnow_Inner(*surf_, *snow_, *met_, *params_, *eb_);
    return eb_->fQm;
  }

 private:
  GroundProperties const* const surf_;
  ModelParams const* const params_;
  MetData const* const met_;

  SnowProperties* const snow_;
  EnergyBalance* const eb_;
};


// Convergence criteria for root-finding
struct Tol_ {
  Tol_(double eps) : eps_(eps) {}
  bool operator()(const double& a, const double& b) const { return std::abs(a - b) <= eps_; }
  double eps_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
