/*
  Functions for calculating the snow-surface energy balance.

  Incoming Longwave radation is cacualted in this version, but if data is
  available we could incorporate it with the available met data.

  Atmospheric pressure is often used in snow models, If data is available we
  could incorperate it but for now Pa is held constant at 100 Pa.

  *** Equation for saturated vapor pressure over water is taken from Bolton,
      1980 'Monthly Weather Review'

  *** Equation for saturated vaport pressure over snow is taken from Buck,
      1996 'Buck Research Manual'

  *** See: http://cires.colorado.edu/~voemel/vp.html

*/

#include <iostream>
#include <cmath>

#include "SnowEnergyBalance_VPL.hh"

#ifdef ENABLE_DBC
#include "dbc.hh"
#endif

void SurfaceEnergyBalance_VPL::UpdateIncomingRadiation(LocalData& seb) {
  // Calculate incoming short-wave radiation
  seb.st_energy.fQswIn = (1 - seb.st_energy.albedo_value) * seb.st_energy.QswIn;

  // Calculate incoming long-wave radiation
    double EmissivityAir = std::pow((10*seb.vp_air.actual_vaporpressure),(seb.st_energy.temp_air/2016));
    EmissivityAir = (1 - std::exp(-EmissivityAir));
    EmissivityAir = 1.08 * EmissivityAir;
  seb.st_energy.fQlwIn = EmissivityAir * seb.st_energy.stephB * std::pow(seb.st_energy.temp_air,4);

  // Calculate D_h, D_e
  if (seb.st_energy.ht_snow>0.02){// roughness length for wind swept snow = 0.005 m
  seb.st_energy.Dhe = (std::pow(seb.st_energy.VKc,2) * seb.st_energy.Us
                       / std::pow(std::log(seb.st_energy.Zr / seb.st_energy.Zo), 2));
}else{// roughness lenght for open field = 0.03 m
  seb.st_energy.Dhe = (std::pow(seb.st_energy.VKc,2) * seb.st_energy.Us
                       / std::pow(std::log(seb.st_energy.Zr / 0.03), 2));
}
}

void SurfaceEnergyBalance_VPL::UpdateIncomingRadiationDerivatives(LocalData& seb) {
  // Calculate incoming short-wave radiation
  seb.st_energy.fQswIn = 0.;

  // Calculate incoming long-wave radiation
  seb.st_energy.fQlwIn = 0.;

  // Calculate D_h, D_e
  if (seb.st_energy.ht_snow>0.02){// roughness length for wind swept snow = 0.005 m
  seb.st_energy.Dhe = (std::pow(seb.st_energy.VKc,2) * seb.st_energy.Us
                       / std::pow(std::log(seb.st_energy.Zr / seb.st_energy.Zo), 2));
}else{// roughness lenght for open field = 0.03 m
  seb.st_energy.Dhe = (std::pow(seb.st_energy.VKc,2) * seb.st_energy.Us
                       / std::pow(std::log(seb.st_energy.Zr / 0.03), 2));
}
}


void SurfaceEnergyBalance_VPL::UpdateEFluxesSnow(LocalData& seb, double T) {
  double Sqig;
  if (seb.st_energy.Us == 0.) {
    Sqig = 0.;
  } else {
    double Ri  = seb.st_energy.gZr * (seb.st_energy.temp_air-T)
                  / (seb.st_energy.temp_air*std::pow(seb.st_energy.Us,2));
  if (Ri >=0){
      Sqig = 1 / (1 + 10*Ri);
  }else{
     Sqig = (1 - 10*Ri);
  } 
  }
  // Calculate outgoing long-wave radiation
  seb.st_energy.fQlwOut = -seb.st_energy.SEs*seb.st_energy.stephB*std::pow(T,4);

  // Calculate sensible heat flux
  seb.st_energy.fQh = seb.st_energy.rowaCp*seb.st_energy.Dhe*Sqig*(seb.st_energy.temp_air-T);

  // Update vapor pressure of snow
  seb.vp_snow.temp = T;
  UpdateVaporPressure(seb.vp_snow);

  // Calculate latent heat flux
  seb.st_energy.fQe = seb.st_energy.rowaLs*seb.st_energy.Dhe*Sqig*0.622
      * (seb.vp_air.actual_vaporpressure-seb.vp_snow.saturated_vaporpressure) / seb.st_energy.Apa;

  // Calculate heat conducted to ground
  double Ks = 2.9e-6 * std::pow(seb.st_energy.density_snow,2);
  seb.st_energy.fQc = Ks * (T-seb.st_energy.temp_ground) / seb.st_energy.ht_snow;
}


// Determine energy available for melting.
double SurfaceEnergyBalance_VPL::CalcMeltEnergy(LocalData& seb) {
  // Melt energy is the balance
  return seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut
      + seb.st_energy.fQh + seb.st_energy.fQe - seb.st_energy.fQc;
}


// Energy balance for no-snow case.
void SurfaceEnergyBalance_VPL::UpdateGroundEnergy(LocalData& seb) {
  seb.st_energy.fQlwOut = -seb.st_energy.SEtun * seb.st_energy.stephB * std::pow(seb.st_energy.temp_ground,4);
  double porosity=1;
  double porrowaLe;
  double Pc=0.0;
  double Pvl = 0.0;
  double Sqig;
  if (seb.st_energy.Us == 0.) {
    Sqig = 0.;
  } else {
    double Ri = seb.st_energy.gZr * (seb.st_energy.temp_air-seb.st_energy.temp_ground)
        / (seb.st_energy.temp_air*std::pow(seb.st_energy.Us,2));
    if (Ri < 0) { // Unstable condition
      Sqig = (1-10*Ri);
    } else { // Stable Condition
      Sqig = (1/(1+10*Ri));
    }
  }

// }else{
  Pc = 101325 - seb.st_energy.surface_pressure;
   
    
    
// }

  seb.st_energy.fQh = seb.st_energy.rowaCp * seb.st_energy.Dhe * Sqig * (seb.st_energy.temp_air - seb.st_energy.temp_ground);

  double Rair = 1/seb.st_energy.Dhe;

  if (seb.st_energy.water_depth > 0.0) {
    // Checking for standing water
    UpdateVaporPressure(seb.vp_ground);
    // Porosity Smoothing function
    if (seb.st_energy.water_depth < 0.02) {
       double PorosityFactor = seb.st_energy.water_depth/0.02;
       porrowaLe = ((seb.st_energy.surface_porosity * (1-PorosityFactor)) + (1*PorosityFactor)) *  seb.st_energy.rowaLe;
       porosity = ((seb.st_energy.surface_porosity * (1-PorosityFactor)) + (1*PorosityFactor));
    } else {// No porosity 
    porrowaLe = seb.st_energy.rowaLe;
    } 

    if (seb.st_energy.water_depth < 0.01) { // smoothing ponded water transition 
    // Average ponded water Latent heat [ Qe ] 
    double Ponded_Qe;
    double Dry_Qe;
    double Qe_Factor = seb.st_energy.water_depth/0.01;
  
    Ponded_Qe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
       * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;
    
    //Calculate vapor pressure lowering for Dry_Qe
    //Pc = 101325 - seb.st_energy.stored_surface_pressure;
    Pvl = seb.vp_ground.saturated_vaporpressure*std::exp(-Pc/(seb.st_energy.density_w*461.52*seb.st_energy.temp_ground));
    // freeze smoothing for Latent heat
    if (seb.st_energy.temp_ground<275){
        if(seb.st_energy.temp_ground>273.0){
          double QLsmooth = -(273.0 - seb.st_energy.temp_ground);
          QLsmooth = QLsmooth/2;

          Dry_Qe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;

          double QLraw = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;
        //  Dry_Qe = QLsmooth * QLraw;

          // Vapor Pressure Lowering
          QLraw = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;
          //Dry_Qe = QLsmooth * QLraw;
          Dry_Qe = QLraw;
       }else{ // Too cold for evaporation!
     //    Dry_Qe=0;
     //    std::cout<<"Too cold for evaporation!  --> BUT STILL EVAPORATING"<<std::endl;
            Dry_Qe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;
       }
    }else{
         Dry_Qe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
          * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;

        //VP Lowering
         Dry_Qe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
          * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;
    }
    
    double GroundFlux = 0.0;
    GroundFlux = seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut
      + seb.st_energy.fQh + Dry_Qe;
    GroundFlux = seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut
      + seb.st_energy.fQh + Ponded_Qe;
    Dry_Qe = (1.-Qe_Factor) * Dry_Qe;
     Ponded_Qe =  Ponded_Qe * Qe_Factor;
    seb.st_energy.fQe = Ponded_Qe + Dry_Qe;
    
    // Assinging Vapor Mass Flux While subtracting the mass of water added to melt/precipitation in ponded condition 
    // ****** IS THIS APPROPRATE IN THE BARE GROUND TO PONDED TRANSITIOINS ?????  ****
    seb.st_energy.SurfaceVaporFlux = seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); // [m/s]     
    seb.st_energy.Mr -= seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); //Accounting for later addition of water to surface source 

    }else{

     seb.st_energy.fQe =  porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
       * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;
    }    

  } else { // no standing water
   UpdateVaporPressure(seb.vp_ground);
   porosity = seb.st_energy.surface_porosity;
   porrowaLe = porosity * seb.st_energy.rowaLe;
   Pvl = seb.vp_ground.saturated_vaporpressure*std::exp(-Pc/(seb.st_energy.density_w*461.52*seb.st_energy.temp_ground));
   double Vaper_direction = seb.vp_air.actual_vaporpressure-Pvl;  //
   // Equation for reduced vapor diffusivity See Sakagucki and Zeng 2009 eqaution (9) and Moldrup et al., 2004. 
   double Clab_Horn_b = 1;
   double Surface_Vap_Diffusion = std::pow((1-(0.0556/porosity)),(2+3*Clab_Horn_b));
   Surface_Vap_Diffusion = 0.000022 * (std::pow(porosity,2)) * Surface_Vap_Diffusion;
   // Sakagucki and Zeng 2009 eqaution (10)
   double cell_dimension = 0.01/2; // This is from cell center to the boundary.
   double VWC = seb.st_energy.saturation_liquid * porosity;
   double L_Rsoil = std::exp(std::pow((1-(VWC/porosity)),5));
   L_Rsoil = cell_dimension * (L_Rsoil -1) * (1/(2.718-1));
   double Rsoil = L_Rsoil/Surface_Vap_Diffusion;
   if(Vaper_direction <= 0){
    }else{
    Rsoil = 0.0;
    }
   // freeze smoothing for Latent heat
     if (seb.st_energy.temp_ground<275){
        if(seb.st_energy.temp_ground>273.0){
          double QLsmooth = -(273.0 - seb.st_energy.temp_ground); 
          QLsmooth = QLsmooth/2;         
 
          seb.st_energy.fQe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;

          double QLraw = porrowaLe * (1/(Rair+Rsoil)) * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;
          //seb.st_energy.fQe = QLsmooth * QLraw;
          seb.st_energy.fQe = QLraw;

          // Vapor Pressure Lowering
          QLraw = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;
          //seb.st_energy.fQe = QLsmooth * QLraw;
          seb.st_energy.fQe = QLraw;

          QLraw = porrowaLe * (1/(Rair+Rsoil)) * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;
          //seb.st_energy.fQe = QLsmooth * QLraw;
          seb.st_energy.fQe = QLraw;

          // Assinging Vapor Mass Flux While subtracting the mass of water added to melt/precipitation in ponded condition
          seb.st_energy.SurfaceVaporFlux = seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); // [m/s]     
          seb.st_energy.Mr -= seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); //Accounting for later addition of water to surface source

        }else{ // Too cold for evaporation!
         // seb.st_energy.fQe=0;
         seb.st_energy.fQe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;
         seb.st_energy.fQe = porrowaLe * (1/(Rair+Rsoil)) * Sqig * 0.622
            * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;
                  // Assinging Vapor Mass Flux While subtracting the mass of water added to melt/precipitation in ponded condition
          seb.st_energy.SurfaceVaporFlux = seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); // [m/s]     
          seb.st_energy.Mr -= seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); //Accounting for later addition of water to surface source
        }
     }else{
      seb.st_energy.fQe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
        * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;

     seb.st_energy.fQe = porrowaLe * (1/(Rair+Rsoil))  * Sqig * 0.622
        * (seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure) / seb.st_energy.Apa;

   
     //VP Lowering
     seb.st_energy.fQe = porrowaLe * seb.st_energy.Dhe * Sqig * 0.622
        * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;

     seb.st_energy.fQe = porrowaLe * (1/(Rair+Rsoil))  * Sqig * 0.622
        * (seb.vp_air.actual_vaporpressure-Pvl) / seb.st_energy.Apa;


     // Assinging Vapor Mass Flux While subtracting the mass of water added to melt/precipitation in ponded condition
     seb.st_energy.SurfaceVaporFlux = seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); // [m/s] 
     seb.st_energy.Mr -= seb.st_energy.fQe / (seb.st_energy.density_w*seb.st_energy.Le); //Accounting for later addition of water to surface source 
     }
  }
// Heat flux to ground surface is the balance.
  seb.st_energy.fQc = seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut
      + seb.st_energy.fQh + seb.st_energy.fQe;
//}else{ // implicit loop
//}

}


// Energy balance for no-snow case.
void SurfaceEnergyBalance_VPL::UpdateGroundEnergyDerivatives(LocalData& seb) {
  seb.st_energy.fQlwOut = -4 * seb.st_energy.SEtun * seb.st_energy.stephB * std::pow(seb.st_energy.temp_ground,3);

  double Sqig, dSqig;
  if (seb.st_energy.Us == 0.) {
    Sqig = 0.;
    dSqig = 0.;
  } else {
    double Ri = seb.st_energy.gZr * (seb.st_energy.temp_air-seb.st_energy.temp_ground)
        / (seb.st_energy.temp_air*std::pow(seb.st_energy.Us,2));
    double dRi = -seb.st_energy.gZr / (seb.st_energy.temp_air*std::pow(seb.st_energy.Us,2));
    if (Ri < 0) { // Unstable condition
      Sqig = (1-10*Ri);
      dSqig = -10*dRi;
    } else { // Stable Condition
      Sqig = (1/(1+10*Ri));
      dSqig = -std::pow(1+10*Ri,-2) * 10 * dRi;
    }
  }

  seb.st_energy.fQh = - seb.st_energy.rowaCp * seb.st_energy.Dhe * Sqig
      + seb.st_energy.rowaCp * seb.st_energy.Dhe * dSqig * (seb.st_energy.temp_air - seb.st_energy.temp_ground);

  // -- SKIPPING THIS TERM!
  if (seb.st_energy.water_depth > 0.0) {
    // Checking for standing water
    UpdateVaporPressure(seb.vp_ground);
    seb.st_energy.fQe = 0.;
  } else {
    // no standing water
    seb.st_energy.fQe = 0.;
  }

  // Heat flux to ground surface is the balance.
  seb.st_energy.fQc = seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut
      + seb.st_energy.fQh + seb.st_energy.fQe;

  std::cout << "Energy summary:" << std::endl
            << "  dfQswIn  = " << seb.st_energy.fQswIn << std::endl
            << "  dfQlwIn  = " << seb.st_energy.fQlwIn << std::endl
            << "  dfQlwOut = " << seb.st_energy.fQlwOut << std::endl
            << "  dfQh (s) = " << seb.st_energy.fQh << std::endl
            << "  dfQe (l) = " << seb.st_energy.fQe << std::endl
            << "  dfQc (c) = " << seb.st_energy.fQc << std::endl;
}


// Calculate saturated and actual vapor pressures
void SurfaceEnergyBalance_VPL::UpdateVaporPressure(VaporPressure& vp) {
  double temp;
  //Convert from Kelvin to Celsius
  temp = vp.temp-273.15;
  // Sat vap. press o/water Dingman D-7 (Bolton, 1980)
// *** (Bolton, 1980) Calculates vapor pressure in [kPa]  ****
//  vp.saturated_vaporpressure = 0.6112*std::exp(17.67*temp / (temp+243.5));
  vp.saturated_vaporpressure = (temp>0) ? 0.611*std::exp(17.62*temp/(temp+243.12)) : 0.611*std::exp(22.46*temp/(temp+272.62));
  // (Bolton, 1980)
  vp.actual_vaporpressure = vp.saturated_vaporpressure * vp.relative_humidity;
  // Find dewpoint Temp Dingman D-11
  vp.dewpoint_temp = (std::log(vp.actual_vaporpressure) + 0.4926) / (0.0708-0.00421*std::log(vp.actual_vaporpressure));
  // Convert Tdp from Celsius to Kelvin
  vp.dewpoint_temp = vp.dewpoint_temp + 273.15;
}


// Take a weighted average to get the albedo.
double SurfaceEnergyBalance_VPL::CalcAlbedo(EnergyBalance& eb) {
  double perSnow = 0.0, perTundra=0.0, perWater=0.0;
  // Tundra albedo comes from Grenfell and Perovich, (2004)
  // Water albedo from Cogley J.G. (1979)
  // Albedo for deteriorated ice from Grenfell and Perovich, (2004)
  double AlTundra=0.15, AlWater=0.141, Alice=0.44;
  double AlSnow = 0.0;
  double TransitionVal = eb.AlbedoTrans;  // Set to 2 cm
    double TransionPercent=0.0;
//Shortwave Pentration Depth varries greatly depending on water clarity. 
    double QswPenitrationDepth = 0.1;// This number could easily change

    AlWater=(AlWater*eb.water_fraction) + (Alice*(1-eb.water_fraction));

  if (eb.density_snow <= 432.238) {
    AlSnow = 1.0 - 0.247 * std::pow(0.16 + 110*std::pow(eb.density_snow/1000, 4), 0.5);
  } else {
    AlSnow = 0.6 - eb.density_snow / 4600;
  }

  if (eb.ht_snow > TransitionVal) {
    // Snow is too deep for albedo weighted average, just use all snow
    perSnow = 1.;
  } else if (eb.water_depth <= 0.0) {  // dry ground
    // Transition to dry ground
    perSnow = std::pow(eb.ht_snow/TransitionVal, 2);
    perTundra = 1 - perSnow;
  } else {
    // Transitions to surface water
    perSnow = eb.ht_snow / TransitionVal;
    perWater = 1 - perSnow;
      //Transition from Ponded water to Bare Ground
      if (eb.water_depth < QswPenitrationDepth) {
          TransionPercent=eb.water_depth/QswPenitrationDepth;
          AlWater=AlWater*TransionPercent + AlTundra*(1-TransionPercent);
      }
  }

#ifdef ENABLE_DBC
  AMANZI_ASSERT(std::abs((perSnow + perTundra + perWater) - 1.) < 1.e-16);
#endif

  // weighted average function for surface albedo
  return AlSnow*perSnow + AlTundra*perTundra + AlWater*perWater;
}


// Surface Energy Balance residual
double SurfaceEnergyBalance_VPL::EnergyBalanceResidual(LocalData& seb, double Xx) {
  UpdateEFluxesSnow(seb, Xx);

  // energy balance
  double res = seb.st_energy.ht_snow
      * (seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut
         + seb.st_energy.fQh + seb.st_energy.fQe - seb.st_energy.fQc);
  return res;
}


// Use a bisection method to calculate the temperature of the snow.
double SurfaceEnergyBalance_VPL::CalcSnowTemperature(LocalData& seb) {
  double tol = 1.e-6;
  double deltaX = 5;

  double Xx = seb.st_energy.temp_air;
  double FXx = EnergyBalanceResidual(seb, Xx);
  // NOTE: decreasing function
  // Bracket the root by (a,b)
  double a,b,Fa,Fb;
  if (FXx > 0) {
    b = Xx;
    Fb = FXx;
    a = Xx;
    Fa = FXx;
    while (Fa > 0) {
      b = a;
      Fb = Fa;
      a += deltaX;
      Fa = EnergyBalanceResidual(seb,a);
    }
  } else {
    a = Xx;
    Fa = FXx;
    b = Xx;
    Fb = FXx;
    while (Fb < 0) {
      a = b;
      Fa = Fb;
      b -= deltaX;
      Fb = EnergyBalanceResidual(seb,b);
    }
  }

#ifdef ENABLE_DBC
  AMANZI_ASSERT(Fa*Fb < 0);
#endif

  int maxIterations = 200;
  double res;
  int iter;
  // Bisection Iterations Loop: Solve for Ts using Energy balance equation
  for (int i=0; i<maxIterations; ++i) {
    Xx = (a+b)/2;
    res = EnergyBalanceResidual(seb, Xx);

    if (res>0) {
      b=Xx;
    } else {
      a=Xx;
    }
    if (std::abs(res)<tol) {
      break;
    }
    iter=i;
  }

#ifdef ENABLE_DBC
  AMANZI_ASSERT(std::abs(res) <= tol);
#endif
  return Xx;
}


// Alter mass flux due to melting.
void SurfaceEnergyBalance_VPL::UpdateMassMelt(EnergyBalance& eb) {
  // Melt rate given by energy rate available divided by heat of fusion.
  double melt = eb.Qm / (eb.density_w * eb.Hf);
  eb.Mr += melt;
  eb.MIr -= melt;
}


// Alter mass fluxes due to sublimation/condensation.
void SurfaceEnergyBalance_VPL::UpdateMassSublCond(EnergyBalance& eb) {
  double SublR = -eb.fQe / (eb.density_w * eb.Ls); // [m/s]

  if (SublR < 0 && eb.temp_snow == 273.15) {
    // Condensation, not sublimation.
    // Snow is melting, surface temp = 0 C and condensation is applied as
    // water and drains through snow.  Therefore add directly to melt.
    eb.Mr += -SublR;
  } else {
    //    if (SublR > 0) {
      eb.MIr += -SublR;
      //    }
  }
}


// Alter mass fluxes due to evaporation from soil (no snow).
void SurfaceEnergyBalance_VPL::UpdateMassEvap(EnergyBalance& eb) {
  eb.Mr += eb.fQe / (eb.density_w*eb.Le); // [m/s]
}


// Alter mass fluxes in the case of all snow disappearing.
void SurfaceEnergyBalance_VPL::WaterMassCorrection(EnergyBalance& eb) {
  if (eb.MIr < 0) {
    // convert ht_snow to SWE
    double swe = eb.ht_snow * eb.density_snow / eb.density_w;
    double swe_change = (eb.MIr * eb.dt) + eb.Ps;
    if (swe + swe_change < 0) {
      // No more snow!  Take the rest out of the ground.
      // -- AA re-visit: should we take some from sublimation?
      eb.Mr += (swe + swe_change) / eb.dt;
    }
  }
}


// Calculate snow change ~> settling of previously existing snow (Martinec, 1977)
void SurfaceEnergyBalance_VPL::UpdateSnow(EnergyBalance& eb) {
  if (eb.MIr < 0.) {
    // sublimation, remove snow now
    eb.ht_snow = eb.ht_snow + (eb.MIr * eb.dt * eb.density_w / eb.density_snow);
  }

  // settle the pre-existing snow
  eb.age_snow += eb.dt / 86400.;
  if (eb.age_snow<0){
      eb.age_snow=0;
  }
  double ndensity = std::pow(eb.age_snow,0.3);
//  ndensity = 1;                                                                                                        //TAKE-OUT-AA 
  if (ndensity < 1){// Formula only works from snow older the 1 day
     ndensity = 1;
   }
  double dens_settled = eb.density_freshsnow*ndensity;
//  eb.density_snow = 250;                                                                                                //TAKE-OUT-AA
  double ht_settled = eb.ht_snow * eb.density_snow / dens_settled;

  // Match Frost Age with Assinged density
     //Calculating which Day frost density matched snow Defermation fucntion from (Martinec, 1977) 
  double frost_age = pow((eb.density_frost /eb.density_freshsnow),(1/0.3))-1;
  frost_age = frost_age + eb.dt / 86400.; 

  // determine heights of the sources
  double ht_precip = eb.Ps * eb.density_w / eb.density_freshsnow;
  double ht_frost = eb.MIr > 0. ? eb.MIr * eb.dt * eb.density_w / eb.density_frost : 0.;

  eb.ht_snow = ht_precip + ht_frost + ht_settled;

  // Possibly settling resulted in negative snow pack, if the snow was disappearing?
  eb.ht_snow = std::max(eb.ht_snow, 0.);

  // Take the height-weighted average to determine new density
  if (eb.ht_snow > 0.) {
    eb.density_snow = (ht_precip * eb.density_freshsnow + ht_frost * eb.density_frost
                       + ht_settled * dens_settled) / eb.ht_snow;
//    eb.density_snow = 250;                                                                                             //TAKE-OUT-AA 
  } else {
    eb.density_snow = eb.density_freshsnow;
  }

  // Take the mass-weighted average to determine new age
  if (eb.ht_snow > 0.) {
     eb.age_snow = (eb.age_snow * ht_settled * dens_settled
                     + frost_age * ht_frost * eb.density_frost + eb.dt / 86400. * ht_precip * eb.density_freshsnow)
      / (ht_settled * dens_settled + ht_frost * eb.density_frost + ht_precip * eb.density_freshsnow);    
   } else {
    eb.age_snow = 0;
  }
eb.SWE = eb.ht_snow * eb.density_snow / eb.density_w;
}


// Main snow energy balance function.
void SurfaceEnergyBalance_VPL::SnowEnergyBalance(LocalData& seb) {
  // Caculate Vapor pressure and dewpoint temperature from Air
  UpdateVaporPressure(seb.vp_air);

  // Find effective Albedo
  seb.st_energy.albedo_value = CalcAlbedo(seb.st_energy);
  // Update temperature-independent fluxes, the short- and long-wave incoming
  // radiation.
  UpdateIncomingRadiation(seb);

  // Convert snow precipitation length (m) to SWE length
  //  SWE(seb.st_energy);

  // Initialize mass change rates
  seb.st_energy.Mr = seb.st_energy.Pr / seb.st_energy.dt; // precipitation rain
  seb.st_energy.MIr = 0;
  seb.st_energy.SurfaceVaporFlux = 0.;

  if (seb.st_energy.ht_snow > 0) { // If snow
    // Step 1: Energy Balance
    // Calculate the temperature of the snow.
    seb.st_energy.temp_snow = CalcSnowTemperature(seb);

    // Calculate the energy available for melting, Qm
    if (seb.st_energy.temp_snow <= 273.15) { // Snow is not melting
      seb.st_energy.Qm = 0; //  no water leaving snowpack as melt water
    } else {
      seb.st_energy.temp_snow = 273.15; // Set snow temperature to zero
      UpdateEFluxesSnow(seb, seb.st_energy.temp_snow);
      seb.st_energy.Qm = CalcMeltEnergy(seb); // Recaculate energy balance with melting.
    }

    // Step 2: Mass Balance
    // -- melt
    UpdateMassMelt(seb.st_energy);

    // -- sublimation and condensation rates between ice, melt, and air
    UpdateMassSublCond(seb.st_energy);

    // Make sure proper mass of snowpack water gets delivered to AT
    WaterMassCorrection(seb.st_energy);


  } else { // no snow
    // Energy balance
    UpdateGroundEnergy(seb);

    // Mass balance
    UpdateMassEvap(seb.st_energy);
  }

  // Update snow pack, density
  UpdateSnow(seb.st_energy);

  // set water temp
  seb.st_energy.Trw = seb.st_energy.ht_snow > 0. ? 273.15 : seb.st_energy.temp_air;

  //Final water flux check
   if ((seb.st_energy.ht_snow<0.02)&&(seb.st_energy.Mr<0)){
       seb.st_energy.Mr = 0;
   }

  std::cout << "Energy summary:" << std::endl
            << "  fQswIn  = " << seb.st_energy.fQswIn << std::endl
            << "  fQlwIn  = " << seb.st_energy.fQlwIn << std::endl
            << "  fQlwOut = " << seb.st_energy.fQlwOut << std::endl
            << "  fQh (s) = " << seb.st_energy.fQh << std::endl
            << "  fQe (l) = " << seb.st_energy.fQe << std::endl
            << "  fQc (c) = " << seb.st_energy.fQc << std::endl;
  std::cout << "Mass summary:" << std::endl
            << "  Mr      = " << seb.st_energy.Mr << std::endl
            << "  VapFlux = " << seb.st_energy.SurfaceVaporFlux << std::endl;
  std::cout << "  Snow:" << std::endl
            << "    new ht   = " << seb.st_energy.ht_snow << std::endl
            << "    new age  = " << seb.st_energy.age_snow << std::endl
            << "    new dens = " << seb.st_energy.density_snow << std::endl;

}


// Main energy-only function.
void SurfaceEnergyBalance_VPL::UpdateEnergyBalance(LocalData& seb) {
  if (seb.st_energy.ht_snow > 0.) {
    // // Caculate Vapor pressure and dewpoint temperature from Air
    // UpdateVaporPressure(seb.vp_air);

    // // Find effective Albedo
    // seb.st_energy.albedo_value = CalcAlbedo(seb.st_energy);

    // // Update temperature-independent fluxes, the short- and long-wave incoming
    // // radiation.
    // UpdateIncomingRadiation(seb);

    // seb.st_energy.temp_snow = CalcSnowTemperature(seb);

    // if (seb.st_energy.temp_snow <= 273.15) { // Snow is not melting
    //   seb.st_energy.Qm = 0; //  no water leaving snowpack as melt water
    // } else {
    //   seb.st_energy.temp_snow = 273.15; // Set snow temperature to zero
    //   UpdateEFluxesSnow(seb, seb.st_energy.temp_snow);
    // }

    double Ks = 2.9e-6 * std::pow(seb.st_energy.density_snow,2);
    seb.st_energy.fQc = Ks * (seb.st_energy.temp_snow - seb.st_energy.temp_ground) / seb.st_energy.ht_snow;
  } else {
    // Caculate Vapor pressure and dewpoint temperature from Air
    UpdateVaporPressure(seb.vp_air);

    // Find effective Albedo
    seb.st_energy.albedo_value = CalcAlbedo(seb.st_energy);

    // Update temperature-independent fluxes, the short- and long-wave incoming
    // radiation.
    UpdateIncomingRadiation(seb);

    // Energy balance
    UpdateGroundEnergy(seb);
  }
}


// Main energy-only function.
void SurfaceEnergyBalance_VPL::UpdateEnergyBalanceDerivative(LocalData& seb) {
  if (seb.st_energy.ht_snow > 0.) {
    double Ks = 2.9e-6 * std::pow(seb.st_energy.density_snow,2);
    seb.st_energy.fQc = -Ks / seb.st_energy.ht_snow;
  } else {
    // Caculate Vapor pressure and dewpoint temperature from Air
    UpdateVaporPressure(seb.vp_air);

    // Find effective Albedo
    seb.st_energy.albedo_value = CalcAlbedo(seb.st_energy);

    // Update temperature-independent fluxes, the short- and long-wave incoming
    // radiation.
    UpdateIncomingRadiationDerivatives(seb);

    // Energy balance
    UpdateGroundEnergyDerivatives(seb);
  }
}

