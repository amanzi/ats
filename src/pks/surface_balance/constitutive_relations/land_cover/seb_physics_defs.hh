/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Data structures and parameters for calculating the snow / surface energy balance.

  Nomenclature: A few names are used which imply various surface processes.
    -- snow, water, ice, and tundra are the four possible surface materials
    -- pond indicates ponded water, either ice or water
    -- ground indicates whatever is below any snow -- ice or water or tundra
*/

#ifndef SURFACEBALANCE_SEB_PHYSICS_DEFS_HH_
#define SURFACEBALANCE_SEB_PHYSICS_DEFS_HH_

#include "Teuchos_ParameterList.hpp"
#include "seb_nan.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// needed constants
static const double c_stephan_boltzmann = 0.00000005670373; // Stephan-Boltzmann constant
static const double c_von_Karman = 0.41; // [-] Von Karman constant
static const double c_R_ideal_gas = 461.52; // [Pa m^3 kg^-1 K^-1]

// Catch-all of leftover parameters.
//
// Typically these are thought of as constants in this model, but they may
// change, may be set by the user, or otherwise could vary conceptually.
//
// In particular, we set these in other parts of ATS's input spec, so they can
// be overridden, in which case they really _should_ change to be consistent
// with their usage across the code.  For instance, gravity can be set in ATS's
// input file, so it should be overwritten with that value.  And density_water
// is a field in most of ATS, so this value should be overwritten on a
// grid-cell by grid-cell basis.
struct ModelParams {
  ModelParams() :
      density_air(1.275),       // [kg/m^3]
      density_freshsnow(100.),  // [kg/m^3]
      density_frost(200.),      // [kg/m^3]
      density_snow_max(325.),   // [kg/m^3] // based on observations at Barrow, AK
      thermalK_freshsnow(0.029),// thermal conductivity of fresh snow [W/m K]
      thermalK_snow_exp(2),     // exponent in thermal conductivity of snow model [-]
      H_fusion(333500.0),       // Heat of fusion for melting snow -- [J/kg]
      H_sublimation(2834000.),  // Latent heat of sublimation ------- [J/kg]
      H_vaporization(2497848.), // Latent heat of vaporization ------ [J/kg]
      Cp_air(1004.0),           // Specific heat of air @ const pres- [J/K kg]
      Cv_water(4218.636),       // Specific heat of water ----------- [J/K kg]
      P_atm(101325.),           // atmospheric pressure ------------- [Pa]
      gravity(9.807),           // gravity [kg m / s^2]
      evap_transition_width(100.), // transition on evaporation from surface to
                                   // evaporation from subsurface [m],
                                   // THIS IS DEPRECATED
      Clapp_Horn_b(1.)          // Clapp and Hornberger "b" [-]
  {}

  ModelParams(Teuchos::ParameterList& plist) :
      ModelParams() {
    thermalK_freshsnow = plist.get<double>("thermal conductivity of fresh snow [W m^-1 K^-1]", thermalK_freshsnow);
    thermalK_snow_exp = plist.get<double>("thermal conductivity of snow aging exponent [-]", thermalK_snow_exp);
    density_snow_max = plist.get<double>("max density of snow [kg m^-3]", density_snow_max);
    evap_transition_width = plist.get<double>("evaporation transition width [Pa]", evap_transition_width);
  }

  // likely constants
  double gravity;
  double P_atm;
  double density_air;
  double density_freshsnow;
  double density_frost;
  double density_snow_max;
  double thermalK_freshsnow;
  double thermalK_snow_exp;
  double H_fusion, H_sublimation, H_vaporization;
  double Cp_air, Cv_water;
  double Clapp_Horn_b;

  // other parameters
  double evap_transition_width;
};


// Struct of skin data
struct GroundProperties {
  double temp;                          // temperature [K]
  double pressure;                      // [Pa]
  double ponded_depth;                  // [m]
  double porosity;                      // [-]
  double density_w;                     // density [kg/m^3]
  double dz;                            // [m]
  double albedo;                        // [-]
  double emissivity;                    // [-]
  double saturation_gas;                // [-]
  double roughness;                     // [m] surface roughness of a bare domain
  double snow_death_rate;               // [kg/m^2/s] snow that must die this timestep, make it melt!
  double unfrozen_fraction;             // [-] fraction of ground water that is unfrozen
  double water_transition_depth;        // [m] microtopographic relief, smoothing factor between water and bare ground

  GroundProperties() :
      temp(NaN),
      pressure(NaN),
      porosity(NaN),
      density_w(NaN),
      dz(NaN),
      albedo(NaN),
      emissivity(NaN),
      saturation_gas(NaN),
      roughness(NaN),
      snow_death_rate(0.),
      unfrozen_fraction(0.),
      water_transition_depth(0.02)
  {}

  void UpdateVaporPressure();
};


// Struct of snow state
struct SnowProperties {
  double height;                // snow depth [m] (NOT SWE!)
  double density;               // snow density [ kg / m^3 ]
  double temp;                  // snow temperature [K]
  double albedo;                // [-]
  double emissivity;            // [-]
  double roughness;             // [m] surface roughness of a snow-covered domain

  SnowProperties() :
      height(NaN),
      density(NaN),
      temp(NaN),
      albedo(NaN),
      emissivity(NaN),
      roughness(NaN)
  {}
};


// struct of input MetData.
struct MetData {
  double Us;                    // wind speed, [m/s]
  double Z_Us;
  double QswIn;                 // incoming short-wave radiation, [W/m^2]
  double QlwIn;                 // incoming longwave radiaton, [W/m^2]
  double Ps;                    // precip snow, [m (SWE)/s]
  double Pr;                    // precip rain, [m/s]
  double air_temp;              // air temperature [K]
  double relative_humidity;     // relative humidity [-]

  MetData() :
      Us(NaN),
      Z_Us(NaN),
      QswIn(NaN),
      QlwIn(NaN),
      Ps(NaN),
      Pr(NaN),
      air_temp(NaN),
      relative_humidity(NaN) {}
};


// Struct collecting energy balance terms.
struct EnergyBalance {
  // all are [J/ (m^2 s)]
  double fQswIn;        // incoming short-wave radiation
  double fQlwIn;        // incoming long-wave radiation
  double fQlwOut;       // outgoing long-wave radiation
  double fQh;           // sensible heat
  double fQe;           // latent heat
  double fQc;           // heat conducted to ground surface
  double fQm;           // energy available for melting snow
  double error;         // imbalance!

  EnergyBalance() :
      fQswIn(NaN),
      fQlwIn(NaN),
      fQlwOut(NaN),
      fQh(NaN),
      fQe(NaN),
      fQc(NaN),
      fQm(NaN),
      error(NaN)
  {}
};


// Struct collecting mass balance terms.
struct MassBalance {    // all are in [m/s] of WATER, i.e. snow are in SWE
  double Me;    // condensation of water/frost (if positive),
                // sublimation/evaporation of snow/water (if negative)
  double Mm;    // melt rate (positive indicates increasing water, decreasing snow)
  double dt;    // max dt that may be taken to conserve snow swe

  MassBalance() :
      Me(NaN),
      Mm(NaN) {}
};


// Struct collecting final output fluxes
struct FluxBalance {
  double M_surf; // [m/s], mass to surface system
  double E_surf; // [W/m^2], energy to surface system
  double M_subsurf; // [m/s], mass to/from subsurface system
  double E_subsurf; // [W/m^2], energy to/from subsurface system
  double M_snow; // [m/s], mass swe to snow system

  FluxBalance() :
      M_surf(0.),
      E_surf(0.),
      M_subsurf(0.),
      E_subsurf(0.),
      M_snow(0.) {}
};


// Used to calculate surface properties, prior to calling SEB.
struct SurfaceParams {
  double a_tundra, a_water, a_ice;      // albedos
  double e_snow, e_tundra, e_water, e_ice;     // emissivities
  double Zsmooth, Zrough;      // roughness coefs

  SurfaceParams() :
      a_tundra(0.135),           // [-] Grenfell and Perovich, (2004)
      a_water(0.1168),           // [-] Cogley J.G. (1979)
      a_ice(0.44),              // [-] deteriorated ice, Grenfell and Perovich, (2004)
      e_snow(0.98),             // [-] emissivity for snow, From P. ReVelle (Thesis)
      e_tundra(0.92),           // [-] emissivity for tundra, From P. ReVelle
                                //         (Thesis), Ling & Zhang, 2004
      e_water(0.995),           // [-] emissivity of water, EngineeringToolbox.com
      e_ice(0.98),              // [-] emissivity of ice, EngineeringToolbox.com
      Zsmooth(0.005),           // [m]? roughness coef of smooth
      Zrough(0.03) {}           // [m]? roughness coef of rough
};


// Partitioning of weights for smoothing surface properties, output data structure.
struct Partition {
  double perSnow;
  double perWater;
  double perIce;
  double perTundra;

  double Interpolate(double snow, double water, double ice, double tundra) const {
    return snow*perSnow + water*perWater + ice*perIce + tundra*perTundra;
  }
};


// Partitioning of weights for smoothing surface properties.
struct Partitioner {
  double water_pen, snow_pen;   // radiation penetration depths

  Partitioner() :
      water_pen(0.1),           // [m] QswIn penetration depth through water
      snow_pen(0.02) {}         // [m] QswIn penetration depth through snow
  Partitioner(double water_pen_, double snow_pen_) :
      water_pen(water_pen_),    // [m] QswIn penetration depth through water
      snow_pen(snow_pen_) {}    // [m] QswIn penetration depth through snow

  Partition CalcPartition(double ht_snow, double ht_pond,
                          double unfrozen_fraction) const;
};


} // namespace
} // namespace
} // namespace


#endif
