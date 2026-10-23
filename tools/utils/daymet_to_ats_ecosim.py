"""Download DayMet data and convert it to hourly ATS forcing for ATS-EcoSIM.

Like daymet_to_ats.py (same download, same dataset names, so an ATS input
only changes the file name), but the forcing is hourly with a daily cycle,
and the file carries the site and calendar the EcoSIM PK needs.

Hourly forcing (DayMet 365-day calendar; time 0 = 00:00 local solar time of
the start date; each value is stamped at the start of its hour and holds the
mean of that hour):

  incoming shortwave radiation  the DayMet daily energy (srad x dayl) spread
                                over the hours in proportion to the sine of
                                the sun elevation at mid-hour, with EcoSIM's
                                own declination and sun-angle formulas
                                (MiniFuncMod get_sun_declin, WthrMod
                                PrepHourlyWeather; solar noon 12 h), so that
                                EcoSIM's clear-sky limit never cuts it.
                                Daily totals are preserved.
  air temperature               EcoSIM's daily-to-hourly curves (DayMod TAVG/AMP,
                                WthrMod DailyWeather): minimum at sunrise,
                                maximum at solar noon + 3 h, using the
                                previous day's tmax and the next day's tmin
  vapor pressure air            DayMet daily vp, capped at saturation at the
                                hourly air temperature (as EcoSIM caps VPK)
  precipitation rain / snow     DayMet daily total spread evenly over the day,
                                rain if the hourly air temperature >= 0 C
  wind speed                    4 m s-1 (DayMet has no wind; as daymet_to_ats.py)

Site and calendar, read by the EcoSIM PK through "file"/"header" parameters
(one-element datasets; the same values are also written as attributes):

  latitude [deg], longitude [deg], solar noon [h] (12, local solar time),
  start year, start day of year [0-364], start hour [h] (0)
"""

import datetime
import logging
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import daymet_to_ats  # noqa: E402  (download, read and date validation)

SOLAR_NOON = 12.0       # h, local solar time (EcoSIM SolarNoonHour)
TWILGT = 0.06976        # EcoSimConst: sine of the sun elevation at twilight
SOLAR_CONSTANT = 1360.0  # W m-2 (EcoSIM SolConst = 4.896 MJ m-2 h-1)
WIND_SPEED = 4.0        # m s-1


def read_daymet(filename):
    """Reads a DayMet single-pixel text file. Finds the column-name line
    ("year,yday,...") instead of skipping a fixed number of header lines:
    the header length has changed between DayMet versions, and
    daymet_to_ats.read_daymet (7 lines) misreads current files."""
    with open(filename) as fid:
        for nskip, line in enumerate(fid):
            if line.startswith('year,'):
                break
        else:
            raise RuntimeError('no "year,..." column-name line in %s' % filename)
    return np.genfromtxt(filename, skip_header=nskip, names=True, delimiter=',')


def sun_declination(yday):
    """EcoSIM get_sun_declin [deg], yday = day of year 1..365"""
    xi = np.where(yday == 366, 365.5, yday).astype(float)
    return np.sin(np.radians((xi + 100.0) * 0.9863)) * (-23.47)


def sun_sine(lat, yday, hour_mid):
    """sine of the sun elevation, as EcoSIM PrepHourlyWeather (not clipped)"""
    decl = np.radians(sun_declination(yday))
    azi = np.sin(np.radians(lat)) * np.sin(decl)
    dec = np.cos(np.radians(lat)) * np.cos(decl)
    return azi + dec * np.cos(np.pi / 12.0 * (SOLAR_NOON - hour_mid))


def day_length(lat, yday):
    """EcoSIM GetDayLength [h]"""
    decl = np.radians(sun_declination(yday))
    ratio = (np.sin(np.radians(lat)) * np.sin(decl)) / (np.cos(np.radians(lat)) * np.cos(decl))
    dyl = 12.0 * (1.0 + 2.0 / np.pi * np.arcsin(np.clip(TWILGT + ratio, -1.0, 1.0)))
    dyl = np.where(ratio >= 1.0 - TWILGT, 24.0, dyl)
    return np.where(ratio <= -1.0 + TWILGT, 0.0, dyl)


def hourly_temperature(tmax, tmin, dayl):
    """EcoSIM DayMod/DailyWeather hourly air temperature [C], shape (ndays, 24)"""
    tmax_prev = np.concatenate([tmax[:1], tmax[:-1]])
    tmin_next = np.concatenate([tmin[1:], tmin[-1:]])
    tavg1, tavg2, tavg3 = (tmax_prev + tmin) / 2, (tmax + tmin) / 2, (tmax + tmin_next) / 2
    amp1, amp2, amp3 = tavg1 - tmin, tavg2 - tmin, tavg3 - tmin_next
    J = np.arange(1, 25, dtype=float)[None, :]
    D = dayl[:, None]
    noon = SOLAR_NOON
    night_len = noon + 9.0 - D / 2.0
    t_early = tavg1[:, None] + amp1[:, None] * np.sin((J + noon - 3.0) * np.pi / night_len + np.pi / 2)
    t_late = tavg3[:, None] + amp3[:, None] * np.sin((J - noon - 3.0) * np.pi / night_len + np.pi / 2)
    t_day = tavg2[:, None] + amp2[:, None] * np.sin((J - (noon - D / 2.0)) * np.pi / (3.0 + D / 2.0) - np.pi / 2)
    return np.where(J < noon - D / 2.0, t_early, np.where(J > noon + 3.0, t_late, t_day))


def saturation_vapor_pressure(t_c):
    """[Pa], Tetens over water"""
    return 610.78 * np.exp(17.27 * t_c / (t_c + 237.3))


def daymet_to_ats_ecosim(dat, lat):
    """DayMet daily records -> dict of hourly ATS datasets"""
    yday = dat['yday'].astype(int)
    nd = len(dat)
    hour_mid = np.arange(24) + 0.5                                     # mid-hour, hour J covers [J-1, J)

    # shortwave: daily energy distributed by the sine of the sun elevation
    energy = dat['srad_Wm2'] * dat['dayl_s']                           # J m-2 d-1
    sine = np.clip(sun_sine(lat, yday[:, None], hour_mid[None, :]), 0.0, None)
    weight_sum = sine.sum(axis=1)
    sw = np.where(weight_sum[:, None] > 0,
                  energy[:, None] * sine / np.where(weight_sum > 0, weight_sum, 1.0)[:, None] / 3600.0, 0.0)
    clear = SOLAR_CONSTANT * sine
    # On a day whose DayMet total exceeds EcoSIM's top-of-atmosphere total
    # (its simplified astronomy), every daylight hour is above the limit; cap to
    # it, as EcoSIM itself would (RADN = MIN(SolConst*sin, RADN)).
    ratio = energy / np.where(weight_sum > 0, SOLAR_CONSTANT * weight_sum * 3600.0, np.inf)
    ncap = int(np.sum(sw > clear * (1 + 1e-12)))
    if ncap:
        logging.warning('%d days have a daily shortwave total above EcoSIM\'s top-of-atmosphere '
                        'total (max ratio %.3f); capped to it, as EcoSIM would', int(np.sum(ratio > 1)),
                        ratio.max())
        sw = np.minimum(sw, clear)

    t_c = hourly_temperature(dat['tmax_deg_c'], dat['tmin_deg_c'], day_length(lat, yday))
    vp = np.minimum(dat['vp_Pa'][:, None] * np.ones((1, 24)), saturation_vapor_pressure(t_c))
    precip = (dat['prcp_mmday'] / 1.e3 / 86400.)[:, None] * np.ones((1, 24))

    out = dict()
    out['time [s]'] = np.arange(nd * 24) * 3600.0
    out['air temperature [K]'] = (273.15 + t_c).ravel()
    out['incoming shortwave radiation [W m^-2]'] = sw.ravel()
    out['vapor pressure air [Pa]'] = vp.ravel()
    out['precipitation rain [m s^-1]'] = np.where(t_c >= 0, precip, 0.0).ravel()
    out['precipitation snow [m SWE s^-1]'] = np.where(t_c < 0, precip, 0.0).ravel()
    out['wind speed [m s^-1]'] = WIND_SPEED * np.ones(nd * 24)

    # checks: daily totals preserved
    sw_daily = sw.sum(axis=1) * 3600.0
    if not ncap:
        assert np.allclose(sw_daily, energy, rtol=1e-10, atol=1e-6), 'shortwave daily totals'
    assert np.allclose(precip.sum(axis=1) * 3600.0, dat['prcp_mmday'] / 1.e3, rtol=1e-12), 'precipitation'
    return out


def site_values(lat, lon, start):
    """site/calendar values read by the EcoSIM PK (start: datetime.date at time 0)"""
    return {
        'latitude [deg]': float(lat),
        'longitude [deg]': float(lon),
        'solar noon [h]': SOLAR_NOON,
        'start year': float(start.year),
        'start day of year [0-364]': float(start.timetuple().tm_yday - 1),
        'start hour [h]': 0.0,
    }


def write_ats_ecosim(dat, site, attrs, filename):
    logging.info('Writing ATS-EcoSIM file: {}'.format(filename))
    with h5py.File(filename, 'w') as fid:
        for key, val in dat.items():
            fid.create_dataset(key, data=val)
        for key, val in site.items():
            fid.create_dataset(key, data=np.array([val]))
            fid.attrs[key] = val
        for key, val in attrs.items():
            fid.attrs[key] = val
        fid.attrs['calendar'] = '365_day (DayMet)'
        fid.attrs['time convention'] = ('time [s] = 0 at 00:00 local solar time of the start date; '
                                        'values are hourly means stamped at the start of the hour')


def get_argument_parser():
    parser = daymet_to_ats.get_argument_parser()
    parser.description = __doc__
    parser.add_argument('-o', '--output', help='Output file name (default daymet_hourly_<lat>_<lon>.h5)')
    return parser


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    args = get_argument_parser().parse_args()
    start, end = daymet_to_ats.validate_start_end(args.start, args.end)

    if args.raw_file is not None:
        raw_file = args.raw_file
    else:
        daymet_to_ats.download_daymet(args.directory, args.lat, args.lon, start, end)
        raw_file = os.path.join(args.directory, 'daymet_raw_%s.dat' % daymet_to_ats.file_id(args.lat, args.lon))
    daymet = read_daymet(raw_file)

    if not args.download_only:
        # the first record defines time 0
        first = datetime.date(int(daymet['year'][0]), 1, 1) + datetime.timedelta(days=int(daymet['yday'][0]) - 1)
        ats = daymet_to_ats_ecosim(daymet, args.lat)
        site = site_values(args.lat, args.lon, first)
        attrs = daymet_to_ats.daymet_attrs(args.lat, args.lon, start, end)
        out = args.output or os.path.join(args.directory,
                                          'daymet_hourly_%s.h5' % daymet_to_ats.file_id(args.lat, args.lon))
        write_ats_ecosim(ats, site, attrs, out)
    sys.exit(0)
