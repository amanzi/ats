/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include <cmath>

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "exceptions.hh"
#include "errors.hh"
#include "incident_shortwave_radiation_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
IncidentShortwaveRadiationModel::IncidentShortwaveRadiationModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
IncidentShortwaveRadiationModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  daily_avg_ = plist.get<bool>("daily averaged", true);

  lat_ = plist.get<double>("latitude [degrees]");
  if (lat_ < -90 || lat_ > 90) {
    Errors::Message msg("IncidentShortwaveRadiationModel: \"domain-averaged latitude [degrees]\" "
                        "not in valid range [-90,90]");
    Exceptions::amanzi_throw(msg);
  }

  doy0_ = plist.get<int>("day of year at time 0 [Julian days]", 0);
  if (doy0_ < 0 || doy0_ > 364) {
    Errors::Message msg("IncidentShortwaveRadiationModel: \"day of year at time 0 [Julian days]\" "
                        "not in valid range [0,364]");
    Exceptions::amanzi_throw(msg);
  }
}


// main method
double
IncidentShortwaveRadiationModel::IncidentShortwaveRadiation(double slope,
                                                            double aspect,
                                                            double qSWin,
                                                            double time) const
{
  double time_days = time / 86400.0;
  double doy = std::fmod((double)doy0_ + time_days, (double)365);

  int doy_i = std::lround(doy);
  if (doy_i == 365) {
    // can round up!
    doy_i = 0;
    doy = doy - 365.0;
  }

  double rad = 0.0;
  if (daily_avg_) {
    double hour = 12;
    // to keep this function smooth, we interpolate between neighboring days
    if (doy_i < doy) {
      double rad_i = Impl::Radiation(slope, aspect, doy_i, hour, lat_, qSWin);
      int doy_ii = doy_i + 1;
      if (doy_ii > 364) doy_ii = 0;
      double rad_ii = Impl::Radiation(slope, aspect, doy_ii, hour, lat_, qSWin);
      rad = rad_i + (doy - doy_i) * rad_ii;
    } else {
      double rad_i = Impl::Radiation(slope, aspect, doy_i, hour, lat_, qSWin);
      int doy_ii = doy_i - 1;
      if (doy_ii < 0) doy_ii = 364;
      double rad_ii = Impl::Radiation(slope, aspect, doy_ii, hour, lat_, qSWin);
      rad = rad_i + (doy_i - doy) * rad_ii;
    }
  } else {
    double hour = 12.0 + 24 * (doy - doy_i);
    rad = Impl::Radiation(slope, aspect, doy_i, hour, lat_, qSWin);
  }

  return rad;
}

double
IncidentShortwaveRadiationModel::DIncidentShortwaveRadiationDSlope(double slope,
                                                                   double aspect,
                                                                   double qSWin,
                                                                   double time) const
{
  AMANZI_ASSERT(false);
  return 0;
}

double
IncidentShortwaveRadiationModel::DIncidentShortwaveRadiationDAspect(double slope,
                                                                    double aspect,
                                                                    double qSWin,
                                                                    double time) const
{
  AMANZI_ASSERT(false);
  return 0;
}

double
IncidentShortwaveRadiationModel::DIncidentShortwaveRadiationDIncomingShortwaveRadiation(
  double slope,
  double aspect,
  double qSWin,
  double time) const
{
  return IncidentShortwaveRadiation(slope, aspect, qSWin, time) / qSWin;
}

namespace Impl {

/*Declination angle

    Parameters
    ----------
    doy : int
      Julian day of the year

    Returns
    -------
    delta : double
      The angle, in radians
    */
double
DeclinationAngle(double doy)
{
  return 23.45 * M_PI / 180.0 * std::cos(2 * M_PI / 365 * (172 - doy));
}

/*Angle of the sun as a function of the time of day

    Parameters
    ----------
    hour : double
       Time, in a 24 hour clock in [0,24)

    Returns
    -------
    tau : angle in radians
    */
double
HourAngle(double hour)
{
  return (hour + 12) * M_PI / 12.0;
}

/*The solar altitude

    Parameters
    ----------
    delta : double
      Declination angle (radians)
    phi : double
      Latitude (radians)
    tau : hour angle (radians)

    Returns
    -------
    alpha : altitude (radians)
    */
double
SolarAltitude(double delta, double phi, double tau)
{
  double alpha =
    std::asin(std::sin(delta) * std::sin(phi) + std::cos(delta) * std::cos(phi) * std::cos(tau));
  if (alpha <= 0.25 * M_PI / 180) {
    // sun is beyond the horizon
    alpha = 0.25 * M_PI / 180;
  }
  return alpha;
}

/*The sun's azhimuth

    Parameters
    ----------
    delta : double
      Declination angle (radians)
    phi : double
      Latitude (radians)
    tau : hour angle (radians)

    Returns
    -------
    phi_sun : altitude (radians)
    */
double
SolarAzhimuth(double delta, double phi, double tau)
{
  double phi_sun =
    std::atan(-std::sin(tau) / (tan(delta) * std::cos(phi) - std::sin(phi) * std::cos(tau)));
  if ((phi_sun >= 0) && (-std::sin(tau) <= 0)) {
    phi_sun += M_PI;
  } else if ((phi_sun <= 0) && (-std::sin(tau) >= 0)) {
    phi_sun += M_PI;
  }
  return phi_sun;
}

/*Geometric reference factor for a flat surface.

    Parameters
    ----------
    alpha : double
      solar altitude (radians)
    phi_sun : double
      sun's azhimuth (radians)

    Returns
    -------
    flat : geometric factor [-]
    */
double
FlatGeometry(double alpha, double phi_sun)
{
  return std::sin(alpha);
}

/*Geometric reference factor of a slope at a given aspect.

    Parameters
    ----------
    slope : double or array_like
      Positive, down-dip slope [radians]
    aspect : double or array_like
      Dip direction, clockwise from N = 0 [radians]
    alpha : double
      solar altitude (radians)
    phi_sun : double
      sun's azhimuth (radians)

    Returns
    -------
    factor : double or array_like

    */
double
SlopeGeometry(double slope, double aspect, double alpha, double phi_sun)
{
  return std::cos(slope) * std::sin(alpha) +
         std::sin(slope) * std::cos(alpha) * std::cos(phi_sun - aspect);
}

/*Returns the geometric factor to multiply times a solar radiation to get a slope-aspect specific value

    Parameters
    ----------
    slope : double or array_like
      Positive, down-dip slope, $-grad z \dot \hat{n}^\perp$,
      where $\hat{n}^\perp$ here refers to the projection of
      the normal onto the x-y plane. [-]
    aspect : double or array_like
      Angle of $\hat{n}^\perp$ in map-view, measured
      clockwise from N = 0 [radians].
    doy : int
      Julian day of the year
    hour : double
      Hour of the day, in 24-hour clock [0,24)
    lat : double
      Latitude [degrees]

    Returns
    -------
    Rslope : double or array_like
      Fraction of the full incident sun on a slope.
    Rflat : double or array_like
      Fraction of the full incident sun on a flat surface.
    */
std::pair<double, double>
GeometricRadiationFactors(double slope, double aspect, int doy, double hour, double lat)
{
#ifdef ASSERT_VALID_INPUT
  AMANZI_ASSERT(365 >= doy);
  AMANZI_ASSERT(doy > 0);
  AMANZI_ASSERT(24 > hour);
  AMANZI_ASSERT(hour >= 0);
  AMANZI_ASSERT(slope >= 0);
  AMANZI_ASSERT(360 > aspect);
  AMANZI_ASSERT(aspect >= 0);
#endif

  double delta = DeclinationAngle(doy);
  double lat_r = M_PI / 180. * lat;
  double tau = HourAngle(hour);

  double alpha = SolarAltitude(delta, lat_r, tau);
  double phi_sun = SolarAzhimuth(delta, lat_r, tau);

  double fac_flat = FlatGeometry(alpha, phi_sun);
  double slope_r = std::atan(slope);
  double fac_slope = SlopeGeometry(slope_r, aspect, alpha, phi_sun);
  return std::make_pair(fac_slope, fac_flat);
}

/*Caculates radiation

    Parameters
    ----------
    slope : double
      Positive, down-dip slope, $-grad z \dot \hat{n}^\perp$,
      where $\hat{n}^\perp$ here refers to the projection of
      the normal onto the x-y plane. [-]
    aspect : double
      Angle of $\hat{n}^\perp$ in map-view, measured
      clockwise from N = 0 [radians].
    doy : int
      Julian day of the year
    hour : double
      Hour of the day, in 24-hour clock [0,24)
    lat : double
      Latitude [degrees]
    qSWin : double
      Incoming radiation on a flat surface. [W/m^2]

    Returns
    -------
    qSWincident : double
      Radiation incident on the cell.  [W/m^2] (or the same as qSWin)
*/
double
Radiation(double slope, double aspect, int doy, double hour, double lat, double qSWin)
{
  auto facs = GeometricRadiationFactors(slope, aspect, doy, hour, lat);
  double fac = facs.first / facs.second;
  if (fac > 6.)
    fac = 6.;
  else if (fac < 0.)
    fac = 0.;
  return qSWin * fac;
}

} //namespace Impl
} //namespace Relations
} //namespace SurfaceBalance
} //namespace Amanzi
