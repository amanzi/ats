/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates shortwave as a function of slope/aspect/etc.
/*!

.. _incident-shortwave-radiation-model-spec:
.. admonition:: incident-shortwave-radiation-model-spec

   * `"latitude [degrees]`" ``[double]`` Latitude of the site.  A single
     typical value is fine for most domains, even relatively large ones
     (e.g. HUC4).
   * `"daily averaged`" ``[bool]`` **true** Compute a daily averaged radiation, as
     opposed to a time-specific, diurnal-cycle-resolving value.
   * `"day of year at time 0 [Julian days]`" ``[int]`` **0** (Jan 1).  ATS has
     no notion of absolute time, so to do things that depend upon planetary
     dynamics we must know what the day of the year is.  Typically this is set
     by your meteorological data -- set this to be equal to the day of year of
     met data's time 0.

*/

#ifndef AMANZI_STREAMLIGHT_MODEL_HH_
#define AMANZI_STREAMLIGHT_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class StreamlightModel {
 public:
  explicit StreamlightModel(Teuchos::ParameterList& plist);
 
  int DaysFromJan1(const std::string& mmdd);

  std::pair<int,double> DoyHourFromSeconds(double time_sec, double days_offset);

  struct SolarAngles {
    double solar_declination;
    double solar_altitude;
    double solar_zenith;
    double solar_azimuth;
  };
  SolarAngles CalcSolarAngles(double doy, double hour, double lat, double lon);

  double CorrectSolarAzimuth(const SolarAngles& solang, double doy, double hour, double lat, double lon);

  double DiffuseCalc(double atm_trns, double regreR, double regreK);

  double IntegFunc(double angle, double d_sza, double lai);

  double DtCalc(double lai);

  struct RadiativeTransfer {
    double rad_diff;
    double rad_beam;
    double rad_diff_PAR;
    double rad_beam_PAR;
    double diff_trans_PAR;
    double beam_trans_PAR;
    double total_trans_PAR;
    double diff_trans;
    double beam_trans;
    double total_trans;
    double total_trans_PAR_ppfd;
  };
  RadiativeTransfer GetRadiativeTransferEstimatesCN1998(const SolarAngles& solang, 
                                                        double doy,
                                                        double qSWin,
                                                        double lai);

  std::pair<double,double> ShadeCalc(double solar_altitude,
                                     double delta,
                                     double ponded_depth,
                                     double water_width);

  std::pair<double,double> GetRiparianStreamShading(double solar_azimuth_corr,
                                                    double ponded_depth,
                                                    double solar_altitude);

  struct EnergyStream {
    double energy_diff_surface_PAR;
    double energy_beam_surface_PAR;
    double energy_total_surface_PAR;
    double energy_diff_surface;
    double energy_beam_surface;
    double energy_total_surface;
    double energy_total_surface_PAR_ppfd; 
    double energy_total_streambed;
    double energy_total_streambed_PAR_ppfd;
  };
  EnergyStream GetEnergyStream(const RadiativeTransfer& radtrans,
                               const SolarAngles& solang,
                               double fraction_shade_bank,
                               double fraction_shade_veg,
                               double qSWin,
                               double ponded_depth);

  struct PARtoGPP {
    double GPP_sw_inc_gO2_m2d;
    double GPP_stream_gO2_m2d;
    double GPP_streambed_gO2_m2d;
  };
  PARtoGPP ConvertPARtoGPP(const EnergyStream& estrm, double qSWin);

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
  bool use_leap_;
  double tz_offset_;
  double x_LAD_;
  double bank_height_;
  double bank_slope_; 
  double bottom_width_; 
  double tree_height_;
  double overhang_height_;
  double overhang_;
  double channel_azimuth_; 
  double turbidity_;
  double absorb_coef_clear_water_; 
  double solar_const_;
  double scat_coef_;
  double light_use_eff_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi

#endif
