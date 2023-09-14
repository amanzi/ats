/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates shortwave as a function of slope/aspect/etc.
/*!

.. _incident_shortwave_radiation_model-spec:
.. admonition:: incident_shortwave_radiation_model-spec

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

#ifndef AMANZI_SURFACEBALANCE_INCIDENT_SHORTWAVE_RADIATION_MODEL_HH_
#define AMANZI_SURFACEBALANCE_INCIDENT_SHORTWAVE_RADIATION_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

namespace Impl {
double
DeclinationAngle(double doy);
double
HourAngle(double hour);
double
SolarAltitude(double delta, double phi, double tau);
double
SolarAzhimuth(double delta, double phi, double tau);
double
FlatGeometry(double alpha, double phi_sun);
double
SlopeGeometry(double slope, double aspect, double alpha, double phi_sun);
std::pair<double, double>
GeometricRadiationFactors(double slope, double aspect, int doy, double hour, double lat);
double
Radiation(double slope, double aspect, int doy, double hr, double lat, double qSWin);
} // namespace Impl


class IncidentShortwaveRadiationModel {
 public:
  explicit IncidentShortwaveRadiationModel(Teuchos::ParameterList& plist);

  double IncidentShortwaveRadiation(double slope, double aspect, double qSWin, double time) const;

  double
  DIncidentShortwaveRadiationDSlope(double slope, double aspect, double qSWin, double time) const;
  double
  DIncidentShortwaveRadiationDAspect(double slope, double aspect, double qSWin, double time) const;
  double DIncidentShortwaveRadiationDIncomingShortwaveRadiation(double slope,
                                                                double aspect,
                                                                double qSWin,
                                                                double time) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:
  bool daily_avg_;
  double lat_;
  int doy0_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi

#endif
