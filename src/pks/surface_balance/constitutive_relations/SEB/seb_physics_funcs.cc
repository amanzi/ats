/*
  Functions for calculating the snow / surface energy balance.
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <boost/math/tools/roots.hpp>

#include "dbc.hh"

#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {

#define SWE_EPS 1.e-12

void UpdateIncomingRadiation(const SEB& seb, EnergyBalance& eb, bool debug, const Teuchos::RCP<VerboseObject>& vo) {
  // Calculate incoming short-wave radiation
  eb.fQswIn = (1 - seb.in.surf.albedo) * seb.in.met.QswIn;

//  // Calculate incoming long-wave radiation
//  const ThermoProperties& vp_air = seb.in.met.vp_air;
//  double e_air = std::pow(10*vp_air.actual_vaporpressure, vp_air.temp / 2016.);
//  e_air = 1.08 * (1 - std::exp(-e_air));
  eb.fQlwIn = seb.in.met.QlwIn;

  // Calculate D_h, D_e, 
  eb.Dhe = std::pow(seb.params.VKc,2) * seb.in.met.Us
                       / std::pow(std::log(seb.params.Zr / seb.in.surf.Zo), 2);
  if (debug) {
    std::cout << "Incoming Radiation Energy Terms:" << "\n"
              << "  windspeed, Zo: " << seb.in.met.Us <<"  "<<seb.in.surf.Zo << "\n"
              << "  fQswIn   = " << eb.fQswIn << "\n"
              << "  fQlwIn   = " << eb.fQlwIn << "\n"
              << "  wind Ref Ht [m] = " << seb.params.Zr << std::endl;
  }    
}

void UpdateEvapResistance(const SEB& seb, const ThermoProperties& vp_surf, EnergyBalance& eb, bool debug, const Teuchos::RCP<VerboseObject>& vo) {
   const ThermoProperties& vp_air = seb.in.met.vp_air;
   
   double Vapor_direction = vp_air.actual_vaporpressure - vp_surf.actual_vaporpressure;  //
   
// Equation for reduced vapor diffusivity See Sakagucki and Zeng 2009 eqaution (9) and Moldrup et al., 2004. 
   double Clab_Horn_b = 1;
   double actual_porosity = 0.9;  // Hard coded for moss Fix this to pass in form ATS ~AA  *******************
   double Surface_Vap_Diffusion = std::pow((1-(0.0556/actual_porosity)),(2+3*Clab_Horn_b));
   Surface_Vap_Diffusion = 0.000022 * (std::pow(actual_porosity,2)) * Surface_Vap_Diffusion;
// Sakagucki and Zeng 2009 eqaution (10)
   double cell_dimension = 0.01/2; // This is from cell center to the boundary **** HARD CODED FOR CURENT MOSS CELLS ********.
   double L_Rsoil = std::exp(std::pow(1-seb.in.surf.saturation_liquid, 5));
   L_Rsoil = cell_dimension * (L_Rsoil -1) * (1/(std::exp(1.)-1));
   double Rsoil = 0.0;
   
   if(Vapor_direction <= 0){
      Rsoil = L_Rsoil/Surface_Vap_Diffusion;
      }else{
      Rsoil = 0.0;
   }
   double Rair = 1/eb.Dhe;
   eb.Evap_Resistance = Rair;

   if((seb.in.snow_old.ht==0)&&(seb.in.surf.saturation_liquid < 1)) {
      eb.Evap_Resistance = Rair + Rsoil;
   }else{  
      eb.Evap_Resistance = Rair;
   }

   if (debug) {
     std::cout<<"Snow ht: "<<seb.out.snow_new.ht<<"  Saturation: "<<seb.in.surf.saturation_liquid<<std::endl;
     std::cout<<"Air_Temp: "<<vp_air.temp<<" Ground Temp: "<<vp_surf.temp<< std::endl;
     std::cout<<"Vapor Direction: "<<Vapor_direction<<" Air Vapor Pres: "<<vp_air.actual_vaporpressure<<" Surf Vapor Pres: "<<vp_surf.actual_vaporpressure << std::endl;
     std::cout<<"Surface_Vap_Diffusion: "<< Surface_Vap_Diffusion <<"  L_Rsoil:  "<<L_Rsoil<<"  Rsoil: "<<Rsoil<<"   sat_l: "<<seb.in.surf.saturation_liquid<<std::endl;
     std::cout<<"Air Resistance: "<<Rair<<std::endl;
     std::cout<<"Evap Resistance: "<<eb.Evap_Resistance<<std::endl; 
   }
   
}


void UpdateEnergyBalance(const SEB& seb, const ThermoProperties& vp_surf, EnergyBalance& eb, bool debug, const Teuchos::RCP<VerboseObject>& vo) {
  const ThermoProperties& vp_air = seb.in.met.vp_air;

  // Calculate outgoing long-wave radiation
  eb.fQlwOut = -seb.in.surf.emissivity*seb.params.stephB*std::pow(vp_surf.temp,4);

  double Sqig;          // Stability function for use in e and h, precalculated for efficiency
  double air_temp = seb.in.met.vp_air.temp;
  double Ri  = seb.params.gravity * seb.params.Zr * (air_temp - vp_surf.temp)
      / (air_temp * std::pow(seb.in.met.Us,2));
  if (Ri >= 0.) {
    // stable condition or snow?
    Sqig = 1 / (1 + 10*Ri);
  } else {
    // Unstable condition
    Sqig = (1-10*Ri);
  }

  // Calculate sensible heat flux
  eb.fQh = seb.params.density_air * seb.params.Cp * eb.Dhe * Sqig * (vp_air.temp - vp_surf.temp);

  // Calculate latent heat flux
  double LatenHeatOf = 0.0; 
  if (seb.in.snow_old.ht > 0.){
    LatenHeatOf = seb.params.Ls;
   }else{
    LatenHeatOf = seb.params.Le;
  } 

  // With New Evaporation resistance term
  eb.fQe = vp_surf.porosity * seb.params.density_air * LatenHeatOf * (1/eb.Evap_Resistance) * Sqig * 0.622
      * (vp_air.actual_vaporpressure - vp_surf.actual_vaporpressure) / seb.params.Apa;
  if (debug) {
    std::cout<<"Porosity: "<<vp_surf.porosity<<"  LatentHeatOf: "<<LatenHeatOf<<"  Sqig: "<<Sqig<<std::endl;
    std::cout<<"Evap_Resistance: "<<eb.Evap_Resistance<<" Air vapor pres: "<<vp_air.actual_vaporpressure<<" Snow vapor pres: "<<vp_surf.actual_vaporpressure<<std::endl;
  }

  // Calculate heat conducted to ground, if snow
  if (seb.in.snow_old.ht > 0.) {
    double Ks = seb.params.snoK * 1e-4 * std::pow(seb.in.snow_old.density,seb.params.snoKX);
    double snow_hoar_density = 0;
    if(seb.in.snow_old.density>150){
      snow_hoar_density = 1/((0.90/seb.in.snow_old.density)+(0.10/150));
      Ks = seb.params.snoK * 1e-4 * std::pow(snow_hoar_density,seb.params.snoKX);
    }
    eb.fQc = Ks * (vp_surf.temp - seb.in.vp_ground.temp) / seb.in.snow_old.ht;
  }


  if (debug) {
    std::cout << "Energy Balance Terms (ht_snow = " << seb.in.snow_old.ht << "):" << "\n"
              << "  SnowSurfaceTemp  = " << vp_surf.temp << "\n"
              << "  fQlwOut  = " << eb.fQlwOut << "\n"
              << "  fQh      = " << eb.fQh << "\n"
              << "  fQe      = " << eb.fQe << "\n"
              << "  fQc      = " << eb.fQc << std::endl;
  }
}


void UpdateMassBalance(const SEB& seb, MassBalance& mb, EnergyBalance& eb, SnowProperties& snow_new, bool debug, const Teuchos::RCP<VerboseObject>& vo) {
  // this dt is the max timestep that may be taken to conserve snow mass
  mb.dt = seb.in.dt;

  // Melt rate given by available energy rate divided by heat of fusion.
  mb.Mm = eb.fQm / (seb.in.vp_ground.density_w * seb.params.Hf);

  // Condensation rate
  // Calculate latent heat flux
  double LatenHeatOf = 0.0; 
  if (seb.in.snow_old.ht > 0.){
    LatenHeatOf = seb.params.Ls;
   } else {
    LatenHeatOf = seb.params.Le;
  } 
  mb.Me = eb.fQe / (seb.in.vp_ground.density_w * LatenHeatOf);

  // Snow balance
  if (seb.in.snow_old.ht > 0.) {
    double swe_old = seb.in.snow_old.ht * seb.in.snow_old.density / seb.in.vp_ground.density_w;
    double swe_new = swe_old + (seb.in.met.Ps - mb.Mm + mb.Me)*seb.in.dt;
    // First do a pass to ensure we are not melting or sublimation ALL
    // of the available snow.  If so, adjust dt
    if (swe_new < 0.) {
      // Must adjust... we are melting or sublimating all of the available snow.
      ASSERT(seb.in.met.Ps - mb.Mm + mb.Me < 0.);
      mb.dt = -swe_old / (seb.in.met.Ps - mb.Mm + mb.Me);
      swe_new = 0.;
    }

    // All future calculations work with the new dt, which will
    // exactly result in 0 snow height (if it would have been negative).

    // age the old snow
    double age_settled = seb.in.snow_old.age + mb.dt / 86400.;
    double dens_settled = seb.params.density_freshsnow
        * std::max(std::pow(age_settled, 0.3), 1.);

    // Match frost age with assigned density -- Calculate which day frost
    // density matched snow defermation function from (Martinec, 1977)
    double age_frost = std::pow((seb.params.density_frost / seb.params.density_freshsnow),
				(1/0.3)) - 1 + mb.dt / 86400.;

    // precip
    double age_precip = mb.dt / 86400.;

    // determine the new height
    // -- sources
    double swe_settled = swe_old;
    double swe_frost = mb.Me > 0. ? mb.Me*mb.dt : 0.;
    double swe_precip = seb.in.met.Ps*mb.dt;

    // -- sinks
    double swe_subl =  mb.Me < 0. ? -mb.Me*mb.dt : 0.;
    double swe_melt = mb.Mm * mb.dt;

    // -- sublimate precip first
    ASSERT(swe_subl >= 0.);
    if (swe_subl > 0.) {
      if (swe_subl > swe_precip) {
	swe_subl -= swe_precip;
	swe_precip = 0.;
      } else {
	swe_precip -= swe_subl;
	swe_subl = 0.;
      }
    }

    // -- next sublimate settled snow
    if (swe_subl > 0.) {
      ASSERT(swe_subl <= swe_settled + SWE_EPS);
      swe_settled -= swe_subl;
      swe_subl = 0.;
    }

    // -- melt settled snow first
    ASSERT(swe_melt >= -SWE_EPS);
    if (swe_melt > 0.) {
      if (swe_melt > swe_settled) {
        swe_melt -= swe_settled;
        swe_settled = 0.;
      } else {
        swe_settled -= swe_melt;
        swe_melt = 0.;
      }
    }

    // -- now melt frost, precip by even amounts
    if (swe_melt > SWE_EPS) {
      ASSERT(swe_frost + swe_precip > 0.);
      double swe_melt_from_frost = swe_melt * (swe_frost / (swe_frost + swe_precip));
      double swe_melt_from_precip = swe_melt - swe_melt_from_frost;

      swe_frost -= swe_melt_from_frost;
      swe_precip -= swe_melt_from_precip;
    }

    // -- check we didn't screw up
    ASSERT(swe_settled >= -SWE_EPS);
    ASSERT(swe_frost >= -SWE_EPS);
    ASSERT(swe_precip >= -SWE_EPS);

    // -- convert these to heights
    double ht_settled = swe_settled * seb.in.vp_ground.density_w / dens_settled;
    double ht_frost = swe_frost * seb.in.vp_ground.density_w / seb.params.density_frost;
    double ht_precip = swe_precip * seb.in.vp_ground.density_w / seb.params.density_freshsnow;

    // set the snow properties
    double swe_total = std::max(swe_settled + swe_frost + swe_precip, 0.);
    ASSERT(std::abs(swe_total - swe_new) < SWE_EPS);
    snow_new.ht = std::max(ht_settled + ht_frost + ht_precip, 0.);
    snow_new.age = swe_new > 0. ? (swe_settled*age_settled + swe_frost*age_frost + swe_precip*age_precip) / swe_new : 0.;
    snow_new.density = snow_new.ht > 0. ? swe_new * seb.in.vp_ground.density_w / snow_new.ht : seb.params.density_freshsnow;
    snow_new.SWE = swe_new;

    // set the water properties
    // -- water source to ground is (corrected) melt and rainfall
    // NOTE: these rates can only be correct if over mb.dt
    mb.MWg = mb.Mm + seb.in.met.Pr;
    mb.MWg_subsurf = 0.;
    mb.MWg_temp = (mb.MWg > 0. && mb.Mm > 0.) ? (mb.Mm * 273.15 + seb.in.met.Pr * seb.in.met.vp_air.temp) / mb.MWg : (seb.in.met.vp_air.temp > 273.15) ? seb.in.met.vp_air.temp: 273.15;

  } else {
    // set the snow properties
    snow_new.ht = seb.in.met.Ps * seb.in.dt
        * seb.in.vp_ground.density_w / seb.params.density_freshsnow;
    snow_new.age = seb.in.dt / 86400.;
    snow_new.density = seb.params.density_freshsnow;
    snow_new.SWE = snow_new.ht * snow_new.density / seb.in.vp_ground.density_w;

    // set the water properties
    // -- water source to ground is rainfall + condensation
    // -- evaporation is taken from ground if ponded water, from cell source if not (with transition)
    mb.MWg_temp = (seb.in.met.vp_air.temp > 273.15) ? seb.in.met.vp_air.temp: 273.15;
    mb.MWg = seb.in.met.Pr;
    mb.MWg_subsurf = 0.;


    mb.MWg += mb.Me;    
    // if (mb.Me > 0.) {
    //   mb.MWg += mb.Me;

    // } else {
    //   double surf_p = seb.in.vp_ground.pressure;
    //   double trans_factor = surf_p > seb.params.Apa*1000. ? 0. :
    //       surf_p < seb.params.Apa*1000. - seb.params.evap_transition_width ? 1. :
    //       (seb.params.Apa*1000. - surf_p) / seb.params.evap_transition_width;

    //      mb.MWg += (1-trans_factor) * mb.Me;
    //      mb.MWg_subsurf += trans_factor * mb.Me;
    // }
  }

  if (debug) {
    std::cout << "Mass Balance:\n"
              << "  Mm   = " << mb.Mm << "\n"
              << "  Me   = " << mb.Me << "\n"
              << "  Ps   = " << seb.in.met.Ps << "\n"
              << "  Pr   = " << seb.in.met.Pr << "\n"
              << "  Snow Melt:\n"
              << "    old ht   = " << seb.in.snow_old.ht << "\n"
              << "    new ht   = " << snow_new.ht << "\n"
              << "    new age  = " << snow_new.age << "\n"
              << "    new dens = " << snow_new.density << "\n"
              << "    SWE      = " << snow_new.SWE << "\n"
              << "  Water Balance:\n"
              << "    surf src = " << mb.MWg << "\n"
              << "    sub src  = " << mb.MWg_subsurf << std::endl;
  }
}


// Snow temperature calculation.
double DetermineSnowTemperature(const SEB& seb, ThermoProperties& vp_snow,
        EnergyBalance& eb, std::string method) {
  SnowTemperatureFunctor_ func(&seb, &vp_snow, &eb);
  Tol_ tol(1.e-10);
  boost::uintmax_t max_it(50);
  double left, right;
  double res_left, res_right;

  double res_init = func(vp_snow.temp);
  if (res_init < 0.) {
    right = vp_snow.temp;
    res_right = res_init;

    left = vp_snow.temp - 1.;
    res_left = func(left);
    while (res_left < 0.) {
      right = left;
      res_right = res_left;
      left = left - 1.;
      res_left = func(left);
    }
  } else {
    left = vp_snow.temp;
    res_left = res_init;

    right = vp_snow.temp + 1.;
    res_right = func(right);
    while (res_right > 0.) {
      left = right;
      res_left = res_left;
      right = right + 1.;
      res_right = func(right);
    }
  }

  std::pair<double,double> result;
  if (method == "bisection") {
    result = boost::math::tools::bisect(func, left, right, tol, max_it);
  } else if (method == "toms") {
    result = boost::math::tools::toms748_solve(func, left, right, res_left, res_right, tol, max_it);
  }

  return (result.first + result.second)/2.;
}

// master driver
void CalculateSurfaceBalance(SEB& seb, bool debug, const Teuchos::RCP<VerboseObject>& vo) {
  // initialize the data
  seb.in.met.vp_air.UpdateVaporPressure();
  seb.in.vp_ground.UpdateVaporPressure();  // This ground vapor pressure will be ignored in the case of snow

  // Energy balance
  UpdateIncomingRadiation(seb, seb.out.eb, debug, vo);
  
  // Evaporation Resistance Term
  UpdateEvapResistance(seb, seb.in.vp_ground, seb.out.eb, debug, vo);  

  if (seb.in.snow_old.ht > 0.) {
    // snow on the ground, solve for snow temperature
    double T_snow = DetermineSnowTemperature(seb, seb.in.vp_snow, seb.out.eb);
    if (T_snow > 273.15) {
      // limit snow temp to 0, then melt with the remaining energy
      seb.in.vp_snow.temp = 273.15;
      seb.in.vp_snow.UpdateVaporPressure();
      UpdateEnergyBalance(seb, seb.in.vp_snow, seb.out.eb, debug, vo);
      seb.out.eb.BalanceViaMelt();
    } else {
      // snow not melting
      seb.in.vp_snow.temp = T_snow;
      seb.in.vp_snow.UpdateVaporPressure();
      UpdateEnergyBalance(seb, seb.in.vp_snow, seb.out.eb, debug, vo);
      seb.out.eb.fQm = 0.;
    }

  } else {
    // no snow on the ground, balance given by conduction
    seb.in.vp_ground.UpdateVaporPressure();
    UpdateEnergyBalance(seb, seb.in.vp_ground, seb.out.eb, debug, vo);
    seb.out.eb.BalanceViaConduction();
  }

  // Mass balance
  UpdateMassBalance(seb, seb.out.mb, seb.out.eb, seb.out.snow_new, debug, vo);

  if (vo.get()) {
    *vo->os() << "---------------------------------------------------------" << std::endl
              << "Surface Energy Balance:" << std::endl
              << "  Incoming Radiation Energy Terms:" << std::endl
              << "    windspeed, Zo: " << seb.in.met.Us << "  " << seb.in.surf.Zo << std::endl
              << "    fQswIn   = " << seb.out.eb.fQswIn << std::endl
              << "    fQlwIn   = " << seb.out.eb.fQlwIn << std::endl
              << "    wind Ref Ht [m] = " << seb.params.Zr << std::endl
              << "  Evap/Cond Terms:" << std::endl
              << "    air temp, skin temp: " << seb.in.met.vp_air.temp << "  " << (seb.in.snow_old.ht > 0 ? seb.in.vp_snow.temp : seb.in.vp_ground.temp) << std::endl
              << "    skin saturation: " << seb.in.surf.saturation_liquid << std::endl
              << "    vapor direction (- = evap): " << seb.in.met.vp_air.actual_vaporpressure - (seb.in.snow_old.ht > 0 ? seb.in.vp_snow.actual_vaporpressure : seb.in.vp_ground.actual_vaporpressure) << std::endl
              << "    air,surf vapor pres: " << seb.in.met.vp_air.actual_vaporpressure << "  " << (seb.in.snow_old.ht > 0 ? seb.in.vp_snow.actual_vaporpressure : seb.in.vp_ground.actual_vaporpressure) << std::endl
              << "    evap resistance: " << seb.out.eb.Evap_Resistance << std::endl 
              << "  Energy Balance Terms (ht_snow = " << seb.in.snow_old.ht << "):" << std::endl
              << "    SnowSurfaceTemp  = " << seb.in.vp_snow.temp << std::endl
              << "    GroundSurfaceTemp  = " << seb.in.vp_ground.temp << std::endl
              << "    fQlwOut  = " << seb.out.eb.fQlwOut << std::endl
              << "    fQh      = " << seb.out.eb.fQh << std::endl
              << "    fQe      = " << seb.out.eb.fQe << std::endl
              << "    fQc      = " << seb.out.eb.fQc << std::endl
              << "  Mass Balance:\n"
              << "    Mm   = " << seb.out.mb.Mm << std::endl
              << "    Me   = " << seb.out.mb.Me << std::endl
              << "    Ps   = " << seb.in.met.Ps << std::endl
              << "    Pr   = " << seb.in.met.Pr << std::endl
              << "    Snow Melt:\n"
              << "      old ht   = " << seb.in.snow_old.ht << std::endl
              << "      new ht   = " << seb.out.snow_new.ht << std::endl
              << "      new age  = " << seb.out.snow_new.age << std::endl
              << "      new dens = " << seb.out.snow_new.density << std::endl
              << "      SWE      = " << seb.out.snow_new.SWE << std::endl
              << "    Water Balance:\n"
              << "      surf src = " << seb.out.mb.MWg << std::endl
              << "      sub src  = " << seb.out.mb.MWg_subsurf << std::endl;
  }
}

double CalcAlbedoSnow(double density_snow) {
  double AlSnow;
  if (density_snow <= 432.23309912785146) {
    AlSnow = 1.0 - 0.247 * std::pow(0.16 + 110*std::pow(density_snow/1000, 4), 0.5);
  } else {
    AlSnow = 0.6 - density_snow / 4600;
  }
  return AlSnow;
}

double CalcRoughnessFactor(double air_temp) {
  double Zsmooth = 0.005;
  double Zrough = 0.04;

  double Zfraction = air_temp < 270. ? 1. : air_temp > 280. ? 0. : -0.1*air_temp + 28;
  return Zsmooth * Zfraction + Zrough * (1. - Zfraction);
}

} // namespace
} // namespace
} // namespace
