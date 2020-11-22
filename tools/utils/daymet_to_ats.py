"""Download DayMet data and convert it to ATS format.

DayMet is downloaded in point mode based on lat-lon, then converted to
hdf5 files that ATS knows how to read.
"""

import requests
import datetime
import logging
import h5py
import sys, os
import numpy as np

def file_id(lat, lon):
    """Returns a lat-lon id for use in filenames"""
    return '{:.4f}_{:.4f}'.format(lat,lon).replace('.','p')


def daymet_rest_url(lat, lon, start, end, vars=None):
    """Generates the DayMet Rest API URL."""
    #
    # NOTE: it is unclear why this was implemented this way -- shouldn't it use
    # the standard requests dictionary-style syntax?  But if it ain't broke
    # don't fix it. --etc
    #    
    daymet_vars = ['dayl', 'prcp', 'srad', 'swe', 'tmax', 'tmin', 'vp']    

    # check variable names are valid
    if vars is None:
        vars = daymet_vars
    elif type(vars) is list:
        for v in vars:
            if v not in daymet_vars:
                raise RuntimeError('Requested DayMet variable "%s" is not a valid variable name.  Valid are "%r"'%(v, vald_vars))
    elif type(vars) is str:
        if vars not in daymet_vars:
            raise RuntimeError('Requested DayMet variable "%s" is not a valid variable name.  Valid are "%r"'%(vars, vald_vars))
        vars = [vars,]

    # generate the URL
    base = 'https://daymet.ornl.gov/single-pixel/api/data?'
    lat_str = 'lat={:.4f}'.format(lat)
    lon_str = 'lon={:.4f}'.format(lon)
    var_str = 'vars='+','.join(vars)
    start_str = 'start={}'.format(start)
    end_str = 'end={}'.format(end)

    return base + '&'.join([lat_str, lon_str, var_str, start_str, end_str])

        
def read_daymet(filename):
    """Reads a text file of the form provided by DayMet"""
    return np.genfromtxt(filename, skip_header=7, names=True, delimiter=',')


def download_daymet(outdir, lat, lon, start, end, vars=None):
    """Calls the DayMet Rest API to get data and save raw data."""

    url = daymet_rest_url(lat, lon, start, end, vars)
    logging.info('Querying: %s'%url)
    resp = requests.get(url)

    if resp.status_code != requests.codes.ok:
        logging.warning('  returned code: %r'%resp.status_code)    
        raise RuntimeError('Failed download on "%s" with error code "%r"'%(url, resp.status_code))
    else:
        logging.info('  returned code: %r'%resp.status_code)

    filename = os.path.join(outdir, 'daymet_raw_%s.dat'%file_id(lat, lon))
    logging.info('  writing to disk: %s'%filename)
    with open(filename, 'w') as fid:
        fid.write(resp.text)

    return read_daymet(filename)


def daymet_to_ats(dat):
    """Accepts a numpy named array of DayMet data and returns a dictionary ATS data."""
    dout = dict()
    logging.info('Converting to ATS')

    mean_air_temp_c = (dat['tmin_deg_c'] + dat['tmax_deg_c'])/2.0
    precip_ms = dat['prcp_mmday'] / 1.e3 / 86400.
    
    # Sat vap. press o/water Dingman D-7 (Bolton, 1980)
    sat_vp_Pa = 611.2 * np.exp(17.67 * mean_air_temp_c / (mean_air_temp_c + 243.5))

    dout['time [s]'] = np.arange(0, len(dat), 1)*86400.
    dout['air temperature [K]'] = 273.15 + mean_air_temp_c
    dout['incoming shortwave radiation [W m^-2]'] = dat['srad_Wm2']
    dout['relative humidity [-]'] = np.minimum(1.0, dat['vp_Pa']/sat_vp_Pa)

    dout['precipitation rain [m s^-1]'] = np.where(mean_air_temp_c >= 0, precip_ms, 0)
    dout['precipitation snow [m SWE s^-1]'] = np.where(mean_air_temp_c < 0, precip_ms, 0)

    dout['wind speed [m s^-1]'] = 4. * np.ones_like(dout['time [s]'])
    return dout


def write_ats(dat, attrs, filename):
    """Accepts a dictionary of ATS data and writes it to HDF5 file."""
    logging.info('Writing ATS file: {}'.format(filename))
    with h5py.File(filename, 'w') as fid:
        for key, val in dat.items():
            fid.create_dataset(key, data=val)

        for key, val in attrs.items():
            fid.attrs[key] = val
    return


def get_argument_parser():
    """Gets an argparse parser for use in main"""
    import argparse
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('lat', type=float,
                        help='Latitude, in decimal form, up to 4 decimal digits')
    parser.add_argument('lon', type=float,
                        help='Longitude, in decimal form, up to 4 decimal digits')

    parser.add_argument('-d', '--directory', default='.',
                        help='Directory in which to place output files.')
    
    def string_to_start_date(s):
        if len(s) == 4:
            return datetime.datetime(int(s),1,1)
        else:
            return datetime.datetime.strptime(s, '%Y-%m-%d')
    parser.add_argument('-s', '--start', type=string_to_start_date,
                        help='Start date, either YYYY or YYYY-MM-DD')

    def string_to_end_date(s):
        if len(s) == 4:
            return datetime.datetime(int(s),12,31)
        else:
            return datetime.datetime.strptime(s, '%Y-%m-%d')
    parser.add_argument('-e', '--end', type=string_to_end_date,
                        help='End date, either YYYY or YYYY-MM-DD')

    parser.add_argument('--download-only', action='store_true',
                        help='Only download raw data.')
    parser.add_argument('--raw-file',
                        help='Do not download, and instead use this file as raw data.')
    return parser


def validate_start_end(start, end):
    """Checks that these are valid dates for use with DayMet"""
    daymet_start = datetime.date(1980, 1, 1)
    daymet_end = datetime.date(2018, 12, 31)

    # check start time
    if start is None:
        start = daymet_start
    else:
        if start < daymet_start:
            raise RuntimeError('DayMet starts at "%r", so cannot request data starting at "%r'%(daymet_start, start))

    # check end time
    if end is None:
        end = daymet_end
    else:
        if end > daymet_end:
            raise RuntimeError('DayMet ends at "%r", so cannot request data ending at "%r'%(daymet_end, end))

    # check end and start
    if end <= start:
        raise RuntimeError('Requested start time %r is after requested end time %r'%(start, end))
    return start, end


def daymet_attrs(lat, lon, start, end):
    # set the wind speed height, which is made up
    attrs = dict()
    attrs['wind speed reference height [m]'] = 2.0
    attrs['DayMet latitude [deg]'] = lat
    attrs['DayMet longitude [deg]'] = lon
    attrs['DayMet start date'] = str(start)
    attrs['DayMet end date'] = str(end)
    return attrs


def daymet_filename(lat, lon):
    return 'daymet_raw_%s.h5'%file_id(lat, lon)
    

if __name__ == '__main__':
    parser = get_argument_parser()
    args = parser.parse_args()

    start, end = validate_start_end(args.start, args.end)
    
    # download or read daymet data
    if args.raw_file is not None:
        daymet = read_daymet(args.raw_file)
    else:
        daymet = download_daymet(args.directory, args.lat, args.lon,
                                 start, end)

    # convert to ats and write
    if not args.download_only:
        ats = daymet_to_ats(daymet)
        attrs = daymet_attrs(args.lat, args.lon, start, end)        
        filename_out = daymet_filename(args.lat, args.lon)        
        write_ats(ats, attrs, filename_out)

    sys.exit(0)
    
