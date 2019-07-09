/*
  Data structures and parameters for calculating the snow / surface energy balance.

  Nomenclature: A few names are used which imply various surface processes.
    -- snow, water, ice, and tundra are the four possible surface materials
    -- pond indicates ponded water, either ice or water
    -- ground indicates whatever is below any snow -- ice or water or tundra
*/


#ifndef SURFACEBALANCE_SEB_PHYSICS_DEFS_HH_
#define SURFACEBALANCE_SEB_PHYSICS_DEFS_HH_

#include <limits>
#include "Teuchos_ParameterList.hpp"
#if 0
#define MY_LOCAL_NAN std::numeric_limits<double>::signaling_NaN()
#else
#define MY_LOCAL_NAN std::numeric_limits<double>::quiet_NaN()
#endif



namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {


// Struct of skin data
struct GroundProperties {
  double temp;                          // temperature [K]
  double pressure;                      // [Pa]
  double porosity;                      // [-]
  double density_w;                     // density [kg/m^3]
  double dz;                            // [m]
  double albedo;                        // [-]
  double emissivity;                    // [-]
  double saturation_gas;                // [-]
  double roughness;                     // [m] surface roughness of a bare domain
  double ponded_depth;                  // [m] thickness of ponded water
  double fractional_area;               // [-] not used by SEB, but useful for later bookkeeping
  double snow_death_rate;               // [kg/m^2/s] snow that must die this timestep, make it melt!
  double unfrozen_fraction;		// [-] fraction of ground water that is unfrozen

  GroundProperties() :
      temp(MY_LOCAL_NAN),
      pressure(MY_LOCAL_NAN),
      porosity(MY_LOCAL_NAN),
      density_w(MY_LOCAL_NAN),
      dz(MY_LOCAL_NAN),
      albedo(MY_LOCAL_NAN),
      emissivity(MY_LOCAL_NAN),
      saturation_gas(MY_LOCAL_NAN),
      roughness(MY_LOCAL_NAN),
      fractional_area(0.),
      snow_death_rate(0.),
      unfrozen_fraction(0.)
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
      height(MY_LOCAL_NAN),
      density(MY_LOCAL_NAN),
      temp(MY_LOCAL_NAN),
      albedo(MY_LOCAL_NAN),
      emissivity(MY_LOCAL_NAN),
      roughness(MY_LOCAL_NAN)
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
      Us(MY_LOCAL_NAN),
      Z_Us(MY_LOCAL_NAN),
      QswIn(MY_LOCAL_NAN),
      QlwIn(MY_LOCAL_NAN),
      Ps(MY_LOCAL_NAN),
      Pr(MY_LOCAL_NAN),
      air_temp(MY_LOCAL_NAN),
      relative_humidity(MY_LOCAL_NAN) {}
};


// Catch-all of leftover parameters.
struct ModelParams {

  ModelParams() :
      density_air(1.275),       // [kg/m^3]
      density_water(1000.),     // [kg/m^3]
      density_freshsnow(100.),  // [kg/m^3]
      density_frost(200.),      // [kg/m^3]
      density_snow_max(325.),   // [kg/m^3] // based on observations at Barrow, AK
      thermalK_freshsnow(0.029),// thermal conductivity of fresh snow [W/m K]
      thermalK_snow_exp(2),     // exponent in thermal conductivity of snow model [-]
      Hf(333500.0),             // Heat of fusion for melting snow -- [J/kg]
      Ls(2834000.0),            // Latent heat of sublimation ------- [J/kg]
      Le(2497848.),             // Latent heat of vaporization ------ [J/kg]
      Cp_air(1004.0),           // Specific heat of air @ const pres- [J/K kg]
      Cv_water(4218.636),       // Specific heat of water ----------- [J/K kg]
      VKc(0.41),                // Von Karman Constant -------------- [-]
      stephB(0.00000005670373), // Stephan-Boltzmann constant ------- [W/m^2 K^4]
      Apa(101.325),             // atmospheric pressure ------------- [kPa]
      water_ground_transition_depth(0.02),
      evap_transition_width(100.), // transition on evaporation from surface to evaporation from subsurface [m]
      gravity(9.807),
      Clapp_Horn_b(1.),          // Clapp and Hornberger "b" [-]
      R_ideal_gas(461.52)       // ideal gas law R? [Pa m^3 kg^-1 K^-1]
  {}         // gravity [kg m / s^2]
  
  ModelParams(Teuchos::ParameterList& plist) :
      ModelParams() {
    thermalK_freshsnow = plist.get<double>("thermal conductivity of fresh snow [W m^-1 K^-1]", thermalK_freshsnow);
    thermalK_snow_exp = plist.get<double>("thermal conductivity of snow aging exponent [-]", thermalK_snow_exp);
    density_snow_max = plist.get<double>("max density of snow [kg m^-3]", density_snow_max);
  }
  
  double density_air;
  double density_water;
  double density_freshsnow;
  double density_frost;
  double density_snow_max;
  double thermalK_freshsnow;
  double thermalK_snow_exp;
  
  double Hf, Ls, Le, Cp_air, Cv_water;
  double R_ideal_gas;

  // constants for energy equations
  double VKc;
  double stephB;
  double Clapp_Horn_b;
  
  // other constants
  double Apa;
  double evap_transition_width;
  double water_ground_transition_depth;
  double gravity;

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
      fQswIn(MY_LOCAL_NAN),
      fQlwIn(MY_LOCAL_NAN),
      fQlwOut(MY_LOCAL_NAN),
      fQh(MY_LOCAL_NAN),
      fQe(MY_LOCAL_NAN),
      fQc(MY_LOCAL_NAN),
      fQm(MY_LOCAL_NAN),
      error(MY_LOCAL_NAN)
  {}
};


// Struct collecting mass balance terms.
struct MassBalance {    // all are in [m/s] of WATER, i.e. snow are in SWE
  double Me;    // condensation of water/frost (if positive),
                // sublimation/evaporation of snow/water (if negative)
  double Mm;    // melt rate (positive indicates increasing water, decreasing snow)
  double dt;    // max dt that may be taken to conserve snow swe

  MassBalance() :
      Me(MY_LOCAL_NAN),
      Mm(MY_LOCAL_NAN) {}
};


// Struct collecting final output fluxes
struct FluxBalance {
  double M_surf; // [m/s], mass to surface system
  double E_surf; // [W/m^2], energy to surface system
  double M_subsurf; // [m/s], mass to/from subsurface system
  double E_subsurf; // [W/m^2], energy to/from subsurface system
  double M_snow; // [m/s], mass swe to snow system
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
      e_water(0.995),           // [-] emissivity of water >> tundra, this is from Wiki
      e_ice(0.98),              // [-] emissivity of ice >> tundra, this is from Wiki
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
