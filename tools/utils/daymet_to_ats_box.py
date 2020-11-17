"""Download DayMet data in a box and convert it to Amanzi-style HDF5.

DayMet is downloaded in box mode based on lat-lon bounds, then converted to
hdf5 files that ATS knows how to read.
"""

import requests
import datetime
import logging
import h5py, netCDF4
import sys, os
import numpy as np

VALID_YEARS = (1980,2018)
VALID_VARIABLES = ['tmin', 'tmax', 'prcp', 'srad', 'vp', 'swe', 'dayl']
URL = "http://thredds.daac.ornl.gov/thredds/ncss/grid/ornldaac/1328/{year}/daymet_v3_{variable}_{year}_na.nc4"


def boundsStr(bounds):
    """Returns a lat-lon id for use in filenames"""
    fmt_str = "_".join(['{:,4f}',]*4)
    return fmt_str.format(*bounds).replace('.','p')

def getFilename(tmp, bounds, year, var):
    """Returns a tmp filename for a single download file."""
    filename = 'daymet_{var}_{year}_{north}x{west}_{south}x{east}.nc'.format(year=year, var=var,
                                                                             north=bounds[3], east=bounds[2],
                                                                             west=bounds[0], south=bounds[1])
    return os.path.join(tmp, filename)

def downloadFile(tmp, bounds, year, var, force=False):
    """Gets file for a single year and single variable.

    Parameters
    ----------
    var : str
      Name the variable, see table in the class documentation.
    year : int
      A year in the valid range (currently 1980-2018)
    bounds : [xmin, ymin, xmax, ymax]
      Collect a file that covers this shape or bounds in lat/lon.

    Returns
    -------
    filename : str
      Path to the data file.
    """
    
    if year > VALID_YEARS[1] or year < VALID_YEARS[0]:
        raise ValueError("DayMet data is available from {} to {} (does not include {})".format(VALID_YEARS[0], VALID_YEARS[1], year))
    if var not in VALID_VARIABLES:
        raise ValueError("DayMet data supports variables: {} (not {})".format(', '.join(VALID_VARIABLES), var))

    # get the target filename
    filename = getFilename(tmp, bounds, year, var)

    if not os.path.exists(filename) or force:
        url_dict = {'year':str(year),
                    'variable':var}
        url = URL.format(**url_dict)
        logging.info("  Downloading: {}".format(url))
        logging.info("      to file: {}".format(filename))

        request_params = [('var', 'lat'),
                          ('var', 'lon'),
                          ('var', var),
                          ('west', str(bounds[0])),
                          ('south', str(bounds[1])),
                          ('east', str(bounds[2])),
                          ('north', str(bounds[3])),
                          ('horizStride', '1'),
                          ('time_start', '{}-01-01T12:00:00Z'.format(year)),
                          ('time_end', '{}-12-31T12:00:00Z'.format(year)),
                          ('timeStride', '1'),
                          ('accept', 'netcdf')
                          ]

        r = requests.get(url,params=request_params)
        r.raise_for_status()

        with open(filename, 'wb') as fid:
            fid.write(r.content)
    else:
        logging.info("  Using existing: {}".format(filename))

    return filename
        
class Date:
    """Struct to store day of year and year."""
    def __init__(self, doy, year):
        self.doy = doy
        self.year = year

    def __repr__(self):
        return '{}-{}'.format(self.doy, self.year)

def numDays(start, end):
    """Time difference -- assumes inclusive end date."""
    return 365 * (end.year + 1 - start.year) - (start.doy-1) - (365-end.doy)

def loadFile(fname, var):
    with netCDF4.Dataset(fname, 'r') as nc:
        x = nc.variables['x'][:] * 1000. # km to m
        y = nc.variables['y'][:] * 1000. # km to m
        time = nc.variables['time'][:]
        assert(len(time) == 365)
        val = nc.variables[var][:]
    return x,y,val

def initData(d, vars, num_days, nx, ny):
    for v in vars:
        d[v] = np.zeros((num_days, nx, ny),'d')

def collectDaymet(tmpdir, bounds, start, end, vars=None, force=False):
    """Calls the DayMet Rest API to get data and save raw data."""

    if vars is None:
        vars = VALID_VARIABLES

    dat = dict()
    d_inited = False

    for year in range(start.year, end.year+1):
        for var in vars:
            fname = downloadFile(tmpdir, bounds, year, var, force)
            x,y,v = loadFile(fname, var)
            if not d_inited:
                initData(dat, vars, numDays(start,end), len(x), len(y))
                d_inited = True

            # stuff v in the right spot
            if year == start.year and year == end.year:
                dat[var][:,:,:] = np.transpose(v[start.doy:end.doy+1,:,:], (0,2,1))
            elif year == start.year:
                dat[var][0:365-start.doy+1,:,:] = np.transpose(v[start.doy-1:,:,:], (0,2,1))
            elif year == end.year:
                dat[var][-end.doy:,:,:] = np.transpose(v[-end.doy:,:,:], (0,2,1))
            else:
                my_start = 365 * (year - start.year) - start.doy + 1
                dat[var][my_start:my_start+365,:,:] = np.transpose(v, (0,2,1))

    return dat, x, y

def daymetToATS(dat):
    """Accepts a numpy named array of DayMet data and returns a dictionary ATS data."""
    dout = dict()
    logging.info('Converting to ATS')

    mean_air_temp_c = (dat['tmin'] + dat['tmax'])/2.0
    precip_ms = dat['prcp'] / 1.e3 / 86000. # mm/day --> m/s
    
    # Sat vap. press o/water Dingman D-7 (Bolton, 1980)
    sat_vp_Pa = 611.2 * np.exp(17.67 * mean_air_temp_c / (mean_air_temp_c + 243.5))

    time = np.arange(0, dat['srad'].shape[0], 1)*86400.
    dout['air temperature [K]'] = 273.15 + mean_air_temp_c # K
    dout['incoming shortwave radiation [W m^-2]'] = dat['srad'] # Wm2
    dout['relative humidity [-]'] = np.minimum(1.0, dat['vp']/sat_vp_Pa) # -

    dout['precipitation rain [m s^-1]'] = np.where(mean_air_temp_c >= 0, precip_ms, 0)
    dout['precipitation snow [m SWE s^-1]'] = np.where(mean_air_temp_c < 0, precip_ms, 0)

    dout['wind speed [m s^-1]'] = 4. * np.ones_like(dout['relative humidity [-]'])
    return time, dout

def getAttrs(bounds, start, end):
    # set the wind speed height, which is made up
    attrs = dict()
    attrs['wind speed reference height [m]'] = 2.0
    attrs['DayMet latitude min [deg]'] = bounds[1]
    attrs['DayMet longitude min [deg]'] = bounds[0]
    attrs['DayMet latitude max [deg]'] = bounds[3]
    attrs['DayMet longitude max [deg]'] = bounds[2]
    attrs['DayMet start date'] = str(start)
    attrs['DayMet end date'] = str(end)
    return attrs    

def writeATS(time, dat, x, y, attrs, filename):
    """Accepts a dictionary of ATS data and writes it to HDF5 file."""
    logging.info('Writing ATS file: {}'.format(filename))
    with h5py.File(filename, 'w') as fid:
        fid.create_dataset('time [s]', data=time)
        assert(len(x.shape) == 1)
        assert(len(y.shape) == 1)

        fid.create_dataset('row coordinate [m]', data=x)        
        fid.create_dataset('col coordinate [m]', data=y)

        for key in dat.keys():
            assert(dat[key].shape[0] == time.shape[0])
            assert(dat[key].shape[1] == x.shape[0])
            assert(dat[key].shape[2] == y.shape[0])
            grp = fid.create_group(key)
            for i in range(len(time)):
                grp.create_dataset(str(i), data=dat[key][i,:,:])

        for key, val in attrs.items():
            fid.attrs[key] = val
    return

def validBounds(bounds):
    return True

def getArgumentParser():
    """Gets an argparse parser for use in main"""
    import argparse
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('bounds', type=float, nargs=4,
                        help='Longitude, latitude bounding box: xmin, ymin, xmax, ymax')

    def stringToDate(s):
        if len(s) == 4:
            return Date(1, int(s))

        doy_year = s.split('-')
        if len(doy_year) != 2 or len(doy_year[1]) != 4:
            raise RuntimeError('Invalid date format: {}, should be DOY-YEAR'.format(s))
        return Date(int(doy_year[0]), int(doy_year[1]))
    
    parser.add_argument('-s', '--start', type=stringToDate, default=Date(1,VALID_YEARS[0]),
                        help='Start date, in the form DOY-YEAR (e.g. 274-2018 for Oct 1, 2018)')
    parser.add_argument('-e', '--end', type=stringToDate, default=Date(365,VALID_YEARS[1]),
                        help='End date, in the form DOY-YEAR (e.g. 274-2018 for Oct 1, 2018)')

    parser.add_argument('--download-only', action='store_true',
                        help='Only download raw data.')
    parser.add_argument('--force-download', action='store_true',
                        help='Re-download all files, even if they already exist.')
    parser.add_argument('-o', '--outfile', type=str,
                        help='Output HDF5 filename.')
    return parser



if __name__ == '__main__':
    parser = getArgumentParser()
    args = parser.parse_args()

    validBounds(args.bounds)
    
    if args.outfile is None:
        args.outfile = './daymet_{}_{}_{}.h5'.format(args.start, args.end, boundsStr(args.bounds))

    tmpdir = args.outfile+".tmp"
    os.makedirs(tmpdir, exist_ok=True)

    raw, x, y = collectDaymet(tmpdir, args.bounds, args.start, args.end, args.force_download)

    if not args.download_only:
        time, dout = daymetToATS(raw)
        writeATS(time, dout, x, y, getAttrs(args.bounds, args.start, args.end), args.outfile)

    sys.exit(0)
    
