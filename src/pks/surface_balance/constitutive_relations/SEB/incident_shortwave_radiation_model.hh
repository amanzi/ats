/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates shortwave as a function of slope/aspect/etc.

/*!

Aspect modified shortwave radiation is determined by a factor which
is multiplied by the 'incoming radiation incident on a flat surface'
to determine the 'incoming radiation incident on a sloping surface of
a given aspect' as a function of latitude, slope, aspect, and Julian
day of the year, and time of day.

Note that some careful checking and experimentation has found that, in
general, the daily average incident radiation times the 12-noon aspect
modifier correlates reasonably well with the daily average of the
product of the hourly incident radiation and the hourly aspect
modifier.  It is notably better than the daily average radiation times
the daily average aspect modifier.

Derived from LandLab code, which is released under the MIT license:
https://github.com/landlab/landlab/blob/master/landlab/components/radiation/radiation.py

.. _incident_shortwave_radiation_model-spec:
.. admonition:: incident_shortwave_radiation_model-spec

    * `"daily averaged`" ``[bool]`` **true** Calculate daily averaged values (used for daily averaged input shortwave radiation).
    * `"latitude [degrees]`" ``[double]`` Domain averaged latitude, in degrees.  Must be in the range [-90,90]
    * `"day of year at time 0 [Julian days]`" ``[int]`` **0** Day of the year that the simulation began.  Defaults to 0, or Jan 1.

*/

#ifndef AMANZI_SURFACEBALANCE_INCIDENT_SHORTWAVE_RADIATION_MODEL_HH_
#define AMANZI_SURFACEBALANCE_INCIDENT_SHORTWAVE_RADIATION_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

namespace Impl {
  double DeclinationAngle(double doy);
  double HourAngle(double hour);
  double SolarAltitude(double delta, double phi, double tau);
  double SolarAzhimuth(double delta, double phi, double tau);
  double FlatGeometry(double alpha, double phi_sun);
  double SlopeGeometry(double slope, double aspect, double alpha, double phi_sun);
  std::pair<double,double> GeometricRadiationFactors(double slope, double aspect, int doy, double hour, double lat);
  double Radiation(double slope, double aspect, int doy, double hr, double lat, double qSWin);
}
  

class IncidentShortwaveRadiationModel {

 public:
  explicit
  IncidentShortwaveRadiationModel(Teuchos::ParameterList& plist);

  double IncidentShortwaveRadiation(double slope, double aspect, double qSWin, double time) const;

  double DIncidentShortwaveRadiationDSlope(double slope, double aspect, double qSWin, double time) const;
  double DIncidentShortwaveRadiationDAspect(double slope, double aspect, double qSWin, double time) const;
  double DIncidentShortwaveRadiationDIncomingShortwaveRadiation(double slope, double aspect, double qSWin, double time) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  bool daily_avg_;
  double lat_;
  int doy0_;

};

} //namespace
} //namespace
} //namespace

#endif
