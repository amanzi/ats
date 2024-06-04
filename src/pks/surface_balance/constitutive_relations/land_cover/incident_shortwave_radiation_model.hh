/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates shortwave as a function of slope/aspect/etc.
/*!


Aspect modified shortwave radiation is determined by a factor which is
multiplied by the 'incoming radiation incident on a flat surface' to determine
the 'incoming radiation incident on a sloping surface of a given aspect' as a
function of slope and aspect, Julian day of the year, and time of day.  The
latitude and Julian day of the year are used to modify this with both time of
day and seasonal changes of the planet.

Note that some careful checking and experimentation has found that, in
general, the daily average incoming radiation times the 12-noon aspect
modifier correlates reasonably well with the daily average of the
product of the hourly incoming radiation and the hourly aspect
modifier.  It is notably better than the daily average radiation times
the daily average aspect modifier.

This implementation is derived from `LandLab code
<https://github.com/landlab/landlab/blob/master/landlab/components/radiation/radiation.py>`_,
which is released under the MIT license.

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

    KEYS:
    - `"slope`" **DOMAIN-slope_magnitude**
    - `"aspect`" **DOMAIN-aspect**
    - `"incoming shortwave radiation`" **DOMAIN-incoming_shortwave_radiation**

*/

#pragma once

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

namespace Functions {

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
KOKKOS_INLINE_FUNCTION
double
declinationAngle(double doy)
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
KOKKOS_INLINE_FUNCTION
double
hourAngle(double hour)
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
KOKKOS_INLINE_FUNCTION
double
solarAltitude(double delta, double phi, double tau)
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
KOKKOS_INLINE_FUNCTION
double
solarAzhimuth(double delta, double phi, double tau)
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
KOKKOS_INLINE_FUNCTION
double
flatGeometry(double alpha, double phi_sun)
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
KOKKOS_INLINE_FUNCTION
double
slopeGeometry(double slope, double aspect, double alpha, double phi_sun)
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
KOKKOS_INLINE_FUNCTION
std::pair<double, double>
geometricRadiationFactors(double slope, double aspect, int doy, double hour, double lat)
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

  double delta = declinationAngle(doy);
  double lat_r = M_PI / 180. * lat;
  double tau = hourAngle(hour);

  double alpha = solarAltitude(delta, lat_r, tau);
  double phi_sun = solarAzhimuth(delta, lat_r, tau);

  double fac_flat = flatGeometry(alpha, phi_sun);
  double slope_r = std::atan(slope);
  double fac_slope = slopeGeometry(slope_r, aspect, alpha, phi_sun);
  return std::make_pair(fac_slope, fac_flat);
}

/*

Caculates radiation

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
KOKKOS_INLINE_FUNCTION
double
radiation(double slope, double aspect, int doy, double hr, double lat, double sw_in)
{
  auto facs = geometricRadiationFactors(slope, aspect, doy, hr, lat);
  double fac = facs.first / facs.second;
  if (fac > 6.)
    fac = 6.;
  else if (fac < 0.)
    fac = 0.;
  return sw_in * fac;
}

} // namespace Functions


template <class cView_type, class View_type>
class IncidentShortwaveRadiationModel {
 public:
  static const int n_dependencies = 3;
  static const bool provides_derivatives = false;
  static const std::string eval_type;

  explicit IncidentShortwaveRadiationModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    my_key_ = { Keys::cleanPListName(*plist), Tag{ plist->get<std::string>("tag") } };
    auto domain = Keys::getDomain(my_key_.first);

    slope_key_ =
      Keys::readKeyTag(*plist, domain, "slope magnitude", "slope_magnitude", my_key_.second);
    aspect_key_ = Keys::readKeyTag(*plist, domain, "aspect", "aspect", my_key_.second);
    sw_in_key_ = Keys::readKeyTag(*plist,
                                  domain,
                                  "incoming shortwave radiation",
                                  "incoming_shortwave_radiation",
                                  my_key_.second);

    Teuchos::ParameterList& model_list = plist->sublist("model parameters");
    daily_avg_ = model_list.get<bool>("daily averaged", true);

    lat_ = model_list.get<double>("latitude [degrees]");
    if (lat_ < -90 || lat_ > 90) {
      Errors::Message msg("IncidentShortwaveRadiationModel: \"domain-averaged latitude [degrees]\" "
                          "not in valid range [-90,90]");
      Exceptions::amanzi_throw(msg);
    }
  }

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    sw = res[0];

    slope = deps[0];
    aspect = deps[1];
    sw_in = deps[2];
  }

  void freeViews()
  {
    sw = View_type();
    slope = cView_type();
    aspect = cView_type();
    sw_in = cView_type();
  }

  KeyTagVector getMyKeys() const
  {
    return {
      my_key_,
    };
  }
  KeyTagVector getDependencies() const { return { slope_key_, aspect_key_, sw_in_key_ }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double rad = 0.0;
    if (daily_avg_) {
      double hour = 12;
      // to keep this function smooth, we interpolate between neighboring days
      if (doy_i < doy) {
        double rad_i =
          Functions::radiation(slope(i, 0), aspect(i, 0), doy_i, hour, lat_, sw_in(i, 0));
        int doy_ii = doy_i + 1;
        if (doy_ii > 364) doy_ii = 0;
        double rad_ii =
          Functions::radiation(slope(i, 0), aspect(i, 0), doy_ii, hour, lat_, sw_in(i, 0));
        rad = rad_i + (doy - doy_i) * rad_ii;
      } else {
        double rad_i =
          Functions::radiation(slope(i, 0), aspect(i, 0), doy_i, hour, lat_, sw_in(i, 0));
        int doy_ii = doy_i - 1;
        if (doy_ii < 0) doy_ii = 364;
        double rad_ii =
          Functions::radiation(slope(i, 0), aspect(i, 0), doy_ii, hour, lat_, sw_in(i, 0));
        rad = rad_i + (doy_i - doy) * rad_ii;
      }
    } else {
      double hour = 12.0 + 24 * (doy - doy_i);
      rad = Functions::radiation(slope(i, 0), aspect(i, 0), doy_i, hour, lat_, sw_in(i, 0));
    }
    sw(i, 0) = rad;
  }

  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { assert(false); }
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const { assert(false); }

 public:
  // time variables -- set by the evaluator rather than computed in the loop
  double doy;
  int doy_i;

 protected:
  KeyTag my_key_;

  KeyTag slope_key_;
  KeyTag aspect_key_;
  KeyTag sw_in_key_;

  bool daily_avg_;
  double lat_;

  View_type sw;
  cView_type slope, aspect, sw_in;
};


template <class cView_type, class View_type>
const std::string IncidentShortwaveRadiationModel<cView_type, View_type>::eval_type = "NOT_USED";


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
