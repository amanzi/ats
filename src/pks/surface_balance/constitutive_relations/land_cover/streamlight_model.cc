/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: 
*/

/*

*/

#include "Teuchos_ParameterList.hpp"
#include "streamlight_model.hh"
#include <cmath>

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
StreamlightModel::StreamlightModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
StreamlightModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  // time related parameters
  use_leap_ = plist.get<bool>("consider leap years", false);
  tz_offset_ = plist.get<double>("time zone offset [-]", 0.);

  // channel properties
  x_LAD_ = plist.get<double>("leaf angle distribution [-]", 1.);
  bank_height_ = plist.get<double>("bank height [m]");
  bank_slope_ = plist.get<double>("bank slope [-]"); 
  bottom_width_ = plist.get<double>("channel width at the water-sediment interface [m]"); 
  tree_height_ = plist.get<double>("tree height [m]");
  overhang_height_ = plist.get<double>("height of the maximum canopy overhang [m]"); // height at max canopy radius
  overhang_ = plist.get<double>("effectively max canopy radius [m]");
  channel_azimuth_ = plist.get<double>("channel azimuth [decimal degrees]"); 
  turbidity_ = plist.get<double>("turbidity of water [m^-1]", 0.);
  absorb_coef_clear_water_ = plist.get<double>("absorption coefficient of clear water [m^-1]", 0.1521645); 

  // other
  solar_const_ = plist.get<double>("solar constant for radiative transfer [W m-2]", 1370.);
  scat_coef_ = plist.get<double>("scattering coefficient [-]" , 0.);
  light_use_eff_ = plist.get<double>("light use efficiency [-]" , 0.019); 

}


// Calculate days offset in a year relative to Jan 1st.
int StreamlightModel::DaysFromJan1(const std::string& mmdd)
{
  static const int days_in_month_leap[12] = {31,29,31,30,31,30,31,31,30,31,30,31};
  static const int days_in_month_noleap[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  const int* days_in_month = use_leap_ ? days_in_month_leap : days_in_month_noleap;
  
  int month = std::stoi(mmdd.substr(0,2));
  int day   = std::stoi(mmdd.substr(3,2));
  int days = 0;
  for (int m = 1; m < month; ++m) {
    days += days_in_month[m-1];
  }
  days += day - 1;
  return days;
}


// Calculate day of a year and hour of a day
std::pair<int,double> 
StreamlightModel::DoyHourFromSeconds(double time_sec, double days_offset)
{
  double total_days = time_sec / 86400 + days_offset;
  int days = static_cast<int>(std::floor(total_days));
  int day_of_year = days % 365 + 1;
  double day_fraction = total_days - days;
  double hour_of_day = std::min(static_cast<int>(day_fraction * 24.0) + 1, 24);
  return std::make_pair(day_of_year, hour_of_day);
}


// Calculate all types of solar angles 
StreamlightModel::SolarAngles
StreamlightModel::CalcSolarAngles(double doy, double hour, double lat, double lon)
{
  SolarAngles solang;
  // Julian date
  int jdate = doy - 1 + hour / 24;
  // Solar declination
  solang.solar_declination = 23.45 * M_PI / 180.0 * std::sin(2 * M_PI * (jdate + 284) / 365.25);
  double b = M_PI / 182 * (jdate - 81);
  double equation_of_time = (9.87 * std::sin(2 * b) - 7.53 * std::cos(b) - 1.5 * std::sin(b)) / 1440;
  double mean_solar_time = jdate + (lon - tz_offset_ * 15) / 361;
  double true_solar_time = mean_solar_time + equation_of_time;
  // Solar altitude, adjustment from the Li (2006) code
  double sin_solar_altitude = std::sin(solang.solar_declination) * std::sin(lat * M_PI / 180.0) - 
    std::cos(solang.solar_declination) * std::cos(lat * M_PI / 180.0) * std::cos(2 * M_PI * true_solar_time);
  solang.solar_altitude = std::asin(sin_solar_altitude);
  // Solar zenith
  solang.solar_zenith = M_PI_2 - solang.solar_altitude;
  // Initial estimate of solar azimuth
  double cos_solar_azimuth = std::cos(solang.solar_declination) * 
    std::sin(2 * M_PI * true_solar_time) / std::cos(solang.solar_altitude);
  solang.solar_azimuth = std::acos(std::clamp(cos_solar_azimuth, -1.0, 1.0));
  return solang;
}


// Correct solar azimuth for locations with latitude > solar declination
double 
StreamlightModel::CorrectSolarAzimuth(const SolarAngles& solang, 
                                      double doy,
                                      double hour,
                                      double lat,
                                      double lon)
{
  // Add a small amount of time (1 minute) and recalculate azimuth. 
  // solar_azimuth_shade2
  if (lat * M_PI / 180.0 > solang.solar_declination) {
    SolarAngles solang_tmp = CalcSolarAngles(doy, hour + 1.0 / 60.0, lat, lon);
                                           
    if (solang_tmp.solar_azimuth > solang.solar_azimuth) {
      return M_PI_2 + solang.solar_azimuth;
    } else {
      return M_PI_2 - solang.solar_azimuth;
    }
  } else {
    return M_PI_2 - solang.solar_azimuth;
  }
}


// Calculate diffused fraction, Spitters et al. (1986) Eqs. 20a-20d
// atm_trns: Atmospheric transmission [-]
// R, K: parameters in the regression of diffuse share on transmission [-]
double StreamlightModel::DiffuseCalc(double atm_trns, double regreR, double regreK)
{
  if (atm_trns <= 0.22) {
    return 1.0;
  } else if (atm_trns <= 0.35) {
    return 1.0 - 6.4 * std::pow((atm_trns - 0.22), 2);
  } else if (atm_trns <= regreK) {
    return 1.47 - 1.66 * atm_trns;
  } else {
    return regreR;
  }
}


// Calculate the integral in Eq. 4 of Savoy et al. (2021)
// angle: solar zenith angle [decimal degrees]
// d_sza: differential of solar zenith angle [decimal degrees]
// x_LAD: Leaf angle distribution [-]
double StreamlightModel::IntegFunc(double angle, double d_sza, double lai)
{
  double numerator = std::sqrt(std::pow(x_LAD_, 2) + std::pow(std::tan(angle), 2));
  double denominator = x_LAD_ + 1.774 * std::pow(x_LAD_ + 1.182, -0.733);
  double exponent = std::exp(-numerator / denominator * lai);
  return exponent * std::sin(angle) * std::cos(angle) * d_sza;
}


// Calculate the diffuse transmission coefficient for canopy (C&N (1998) Eq. 15.5)
double StreamlightModel::DtCalc(double lai)
{
  const int N = 90;
  double d_sza = M_PI_2 / N;
  double sum = 0.0;
  for (int deg = 0; deg < N; ++deg) {
    double angle = deg * M_PI / 180.0;
    sum += IntegFunc(angle, d_sza, lai);
  }
  return 2.0 * sum;
}


// Calculate below canopy PAR
// References: Campbell & Norman (1998), Spitters et al. (1986), Goudriaan (1977)
StreamlightModel::RadiativeTransfer
StreamlightModel::GetRadiativeTransferEstimatesCN1998(const SolarAngles& solang, 
                                                      double doy,
                                                      double qSWin,
                                                      double lai)
{
  RadiativeTransfer radtrans;
  // First, partitioning incoming shorwave radiation into beam and diffuse components
  // Following Spitters et al. (1986)
  // ----------------------------------------------------------------------------------
  // Calculate the extra-terrestrial irradiance at a plane parallel 
  // to the earth surface(Spitters et al. (1986) Eq. 1), 
  // i.e., global incoming shortwave radiation [W m-2]
  double Qo = solar_const_ * std::sin(solang.solar_altitude) *
    (1 + 0.033 * std::cos(360 * doy / 365 / 180 * M_PI));

  // The relationship between fraction diffuse and atmospheric transmission from 
  // Spitters et al. (1986) appendix: "The radiation incident upon the earth surface 
  // is partly direct, with angle of incidence equal to the angle of the sun, 
  // and partly diffuse, with incidence under different angles. The diffuse flux 
  // arises from scattering (reflection and transmission) of the sun's rays in the atmosphere. 
  // The share of the diffuse flux will therefore be related to the transmission 
  // of the total radiation through the atmosphere." (Spitters et al., 1986)
  double atm_trns = qSWin / Qo; // Atmospheric transmission
  double regreR = 0.847 - 1.61 * std::sin(solang.solar_altitude) + 
    1.04 * std::pow(std::sin(solang.solar_altitude), 2);
  double regreK = (1.47 - regreR) / 1.66;
  double frac_diff = DiffuseCalc(atm_trns, regreR, regreK);

  // Partition into diffuse and beam radiation
  radtrans.rad_diff = frac_diff * qSWin;           // Diffuse radiation [W m-2]
  radtrans.rad_beam = qSWin - radtrans.rad_diff;   // Beam radiation [W m-2]

  // Second, partition diffuse and beam radiation into PAR following Goudriaan (1977):
  // "Up to now global radiation (300--3000 nm) has been considered, but only 
  // the 400--700nm wavebands are photosynthetically active (PAR); the fraction 
  // PAR amounts to 0.50 and is remarkably constant over different atmospheric 
  // conditions and solar elevation, provided that the angle of sun above 
  // horizon (solar altitude) is > 10 degrees (Szeicz, 1974)." (Spitters et al., 1986)
  // ----------------------------------------------------------------------------------
  radtrans.rad_diff_PAR = 0.5 * radtrans.rad_diff;
  radtrans.rad_beam_PAR = 0.5 * radtrans.rad_beam;

  // Third, calculating beam radiation transmitted through the canopy
  // ----------------------------------------------------------------------------------
  // Calculate the ratio of projected area to hemi-surface area for an ellipsoid 
  // C&N (1998) Eq. 15.4 sensu Campbell (1986)
  double numerator = std::sqrt(std::pow(x_LAD_, 2) + std::pow(std::tan(solang.solar_zenith), 2));
  double denominator = x_LAD_ + 1.774 * std::pow(x_LAD_ + 1.182, -0.733);
  double kbe = numerator / denominator;

  // Fraction of incident beam radiation penetrating the canopy
  // C&N (1998) Eq. 15.1 and leaf absorptivity as 0.8 (C&N (1998) pg. 255) as per Camp
  double tau_b = std::exp(-std::sqrt(0.8) * kbe * lai);

  // Beam radiation transmitted through the canopy
  bool is_night = solang.solar_zenith > M_PI_2;
  radtrans.beam_trans_PAR = is_night ? 0 : radtrans.rad_beam_PAR  * tau_b;
  radtrans.beam_trans = radtrans.beam_trans_PAR / 0.5; 
  
  // Fourth, calculating diffuse radiation transmitted through the canopy
  // ----------------------------------------------------------------------------------
  // Diffuse transmission coefficient for the canopy (C&N (1998) Eq. 15.5)
  double tau_d = DtCalc(lai);
  // Extinction coefficient for black leaves in diffuse radiation
  double kd = -std::log(tau_d) / lai;
  // Diffuse radiation transmitted through the canopy
  radtrans.diff_trans_PAR = is_night ? 0 : radtrans.rad_diff_PAR * std::exp(-std::sqrt(0.8) * kd * lai);
  radtrans.diff_trans = radtrans.diff_trans_PAR / 0.5; 

  // Get the total light transmitted through the canopy
  radtrans.total_trans_PAR = radtrans.diff_trans_PAR + radtrans.beam_trans_PAR; 
  radtrans.total_trans = radtrans.diff_trans + radtrans.beam_trans; 
  radtrans.total_trans_PAR_ppfd = radtrans.total_trans_PAR / 0.235; // convert from W m-2 to umol m-2 s-1

  return radtrans;
}


// Calculate the percent of the wetted width shaded by banks and vegetation
// delta: Difference between the sun and stream azimuth (sun-stream)[decimal degrees]
std::pair<double,double> 
StreamlightModel::ShadeCalc(double solar_altitude,
                            double delta,
                            double ponded_depth,
                            double water_width)
{
  // First, calculating the shading produced by the bank
  // ----------------------------------------------------------------------------------
  // Calculating the length of the shadow perpendicular to the bank produced by the bank
  double bank_shadow_length = (bank_height_ - ponded_depth) * std::sin(delta) / std::tan(solar_altitude);
  // Finding the amount of exposed bank in the horizontal direction
  double exposed_bank = (bank_height_ - ponded_depth) / bank_slope_;
  // Finding how much shade falls on the surface of the water
  double stream_shade_bank = std::max(0., bank_shadow_length - exposed_bank);

  // Second, calculating the shading produced by the Vegetation
  // ----------------------------------------------------------------------------------
  // From top of the tree
  double stream_shade_top = std::max(0., 
    (tree_height_ + bank_height_ - ponded_depth) * std::sin(delta) / 
    std::tan(solar_altitude) - exposed_bank);

  // From the overhang
  double stream_shade_overhang = std::max(0., 
    (overhang_height_ + bank_height_ - ponded_depth) * std::sin(delta) / 
    std::tan(solar_altitude) + overhang_ - exposed_bank);
  
  // Get max(shade from top, shade from overhang)
  double stream_shade_veg_max = std::min(water_width, 
    std::max(stream_shade_top, stream_shade_overhang));
  
  return std::make_pair(stream_shade_bank, stream_shade_veg_max);
}


// Calculate riparian stream shading by SHADE2 model (Li et al. (2012))
std::pair<double,double> 
StreamlightModel::GetRiparianStreamShading(double solar_azimuth_corr,
                                           double ponded_depth,
                                           double solar_altitude)
{
  // First, taking the difference between the sun and stream azimuth (sun-stream)
  // ------------------------------------------------------------------------------
  // This must be handled correctly to determine if the shadow falls towards the river
  // [sin(delta)>0] or towards the bank
  // Eastern shading
  double delta_prime_init = solar_azimuth_corr - channel_azimuth_ * M_PI / 180.;
  double delta_prime = (delta_prime_init >= 0) ? delta_prime_init :
    M_PI + std::fmod(std::abs(delta_prime_init), 2. * M_PI);
  double delta_east = std::fmod(delta_prime, 2. * M_PI);
  // Western shading
  double delta_west = (delta_east >= M_PI) ? delta_east - M_PI : delta_east + M_PI;

  // Second, doing some housekeeping related to bankfull and wetted widths
  // ------------------------------------------------------------------------------
  // Calculate bankfull width
  double bankfull_width = bottom_width_ + bank_height_ / bank_slope_ * 2.;
  // Calculate wetted width
  double water_width = std::min(bankfull_width,
    bottom_width_ + std::min(ponded_depth, bank_height_) / bank_slope_ * 2.);
  
  // Third, calculate the length of shading for each bank
  // ------------------------------------------------------------------------------
  // Calculating shade from the "eastern" bank
  std::pair<double,double> east_bank_veg = ShadeCalc(solar_altitude,
                                                     delta_east, 
                                                     std::min(ponded_depth, bank_height_), 
                                                     water_width);
  double east_bank_shade_length = east_bank_veg.first; 
  double east_veg_shade_length = east_bank_veg.second;
  // Calculating shade from the "western" bank
  std::pair<double,double> west_bank_veg = ShadeCalc(solar_altitude,
                                                     delta_west, 
                                                     std::min(ponded_depth, bank_height_), 
                                                     water_width);
  double west_bank_shade_length = west_bank_veg.first; 
  double west_veg_shade_length = west_bank_veg.second;

  // Fourth, calculate the total length of bank shading
  // ------------------------------------------------------------------------------
  // If total bank shade length is longer than wetted width, set to wetted width
  double total_bank_shade_length = std::min(east_bank_shade_length + west_bank_shade_length, water_width);

  // Fifth, calculate the total length of vegetation shading
  // ------------------------------------------------------------------------------
  // If total vegetation shade length is longer than wetted width, set to wetted width
  double total_veg_shade_length = std::min(east_veg_shade_length + west_veg_shade_length, water_width);

  // Sixth, calculating the percentage of water that is shaded
  // ------------------------------------------------------------------------------
  double fraction_shade_bank = std::min(1., total_bank_shade_length / water_width);
  double fraction_shade_veg = std::min(1., (total_veg_shade_length - total_bank_shade_length) / water_width);

  return std::make_pair(fraction_shade_bank, fraction_shade_veg);
}                                                                                        


// Calculate the weighted mean of light reaching the stream surface
StreamlightModel::EnergyStream
StreamlightModel::GetEnergyStream(const RadiativeTransfer& radtrans,
                                  const SolarAngles& solang,
                                  double fraction_shade_bank,
                                  double fraction_shade_veg,
                                  double qSWin,
                                  double ponded_depth) 
{ 
  EnergyStream estrm;
  bool is_night = solang.solar_zenith > M_PI_2;
  estrm.energy_diff_surface_PAR = is_night ? 0 :
    radtrans.diff_trans_PAR * fraction_shade_veg + 
    radtrans.rad_diff_PAR * (1. - fraction_shade_veg - fraction_shade_bank);

  estrm.energy_beam_surface_PAR = is_night ? 0 : 
    radtrans.beam_trans_PAR * fraction_shade_veg +
    radtrans.rad_beam_PAR * (1. - fraction_shade_veg - fraction_shade_bank);

  estrm.energy_total_surface_PAR = is_night ? 0 :
    radtrans.total_trans_PAR * fraction_shade_veg +
    qSWin * 0.5 * (1. - fraction_shade_veg - fraction_shade_bank);

  estrm.energy_diff_surface = is_night ? 0 :
    radtrans.diff_trans * fraction_shade_veg +
    radtrans.rad_diff * (1. - fraction_shade_veg - fraction_shade_bank);

  estrm.energy_beam_surface = is_night ? 0 :
    radtrans.beam_trans * fraction_shade_veg + 
    radtrans.rad_beam * (1. - fraction_shade_veg - fraction_shade_bank);

  estrm.energy_total_surface = is_night ? 0 :
    radtrans.total_trans * fraction_shade_veg +
    qSWin * (1. - fraction_shade_veg - fraction_shade_bank); 

  estrm.energy_total_surface_PAR_ppfd = is_night ? 0 : 
    radtrans.total_trans_PAR_ppfd * fraction_shade_veg +
    qSWin * 2.114 * (1. - fraction_shade_veg - fraction_shade_bank);

  // Predict the light attenuation coefficient (Kd) based on water turbidity.
  // This method calculates the absorption coefficient (Kd) using an empirical
  // relationship with turbidity. The formula used is derived from a logarithmic 
  // relationship, which estimates the absorption coefficient as a function of
  // turbidity. For cases where turbidity is zero, a default absorption coefficient
  // for clear water is applied. (Savoy, P., & Harvey, J. W. (2021))
  double absorb_coef = (turbidity_ == 0) ? absorb_coef_clear_water_ :
    std::pow(10, (0.52 * std::log10(turbidity_) - 0.26));

  // Calculate the absorption coefficient for clear water and update energy responses.

  // Calculates light absorption as a function of depth in clear water, using the mean 
  // absorption coefficient for pure water derived from 
  // "Pope & Fry 1997 Absorption spectrum ~380–700 nm of pure water. II. Integrating cavity measurements." 
  // Effectively, this assumes that attenuation is solely a function of absorption 
  // from pure water. In other words, the irradiance attenuation coefficient (Kd) 
  // will be equal to just the absorption coefficient for pure water.

  // This method calculates the absorption coefficient (kd) for clear water by
  // adding the absorption coefficient (absorb_coef) and the scattering coefficient (scat_coef). 
  // The scattering coefficient is assumed to be zero due to the lack of a better value, 
  // as per Savoy's assumption. The method then calculates the transmission of light through 
  // the water as a function of depth and updates the energy reaching the streambed.
  double kd = absorb_coef + scat_coef_;
  // Calculate transmission as a function of depth
  double transmission = std::exp(-kd * ponded_depth);
  // Update energy responses for the streambed
  estrm.energy_total_streambed = estrm.energy_total_surface * transmission;
  estrm.energy_total_streambed_PAR_ppfd = estrm.energy_total_surface_PAR_ppfd * transmission;

  return estrm;
}


// Conversion of PAR to mean daily GPP (Savoy, P., & Harvey, J. W. (2021))
StreamlightModel::PARtoGPP
StreamlightModel::ConvertPARtoGPP(const EnergyStream& estrm, double qSWin)
{
  PARtoGPP gpp;
  double common_val = 0.235 * 3600 / 4186.8 * estrm.energy_total_surface * 
    light_use_eff_ * 32 * 0.5 /(12 * 3.4);
  gpp.GPP_sw_inc_gO2_m2d = qSWin / 0.235 * common_val;
  gpp.GPP_stream_gO2_m2d = estrm.energy_total_surface_PAR_ppfd * common_val; 
  gpp.GPP_streambed_gO2_m2d = estrm.energy_total_streambed_PAR_ppfd * common_val;
  return gpp;
}





} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
