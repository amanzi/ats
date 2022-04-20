"""Functionality to use daily meterological data to form smoother daily data, on boxed data."""


import sys, os
import numpy as np
import scipy.signal
import h5py
from matplotlib import pyplot as plt
import logging

_prewritten = ['time [s]', 'x [m]', 'y [m]']

def stack_years(raw):
    """Checks that data can be cast in years and stacks it into a 365 x NYEAR array."""
    ntimes, nx, ny = raw.shape
    if ntimes%365 != 0:
        raise RuntimeError('Raw data provided is not module 365 (len = %d)'%ntimes)

    print(f'read data, nx = {nx}, ny = {ny}, ntimes = {ntimes}')
    nyears = int(np.round(len(raw) / 365))
    raw = raw.reshape((nyears,365,nx,ny))
    return nyears, raw

def smooth_raw_to_spinup(raw, name=None, plot=False):
    """Stacks years and generates a single, smooth year.

    Assumes input met data is daily data, and in 365-day year format such that 
    len(met)%365 == 0.
    """
    nyears, raw = stack_years(raw)
    
    mean = raw.mean(axis=0)
    smooth = scipy.signal.savgol_filter(mean, 101, 3, mode='wrap', axis=0)

    if plot:
        fig, axs = plt.subplots(2,1)
        for i in range(nyears):
            axs[0].plot(raw[i,:,0,0])
        axs[0].set_xticklabels(list())
        axs[1].plot(mean[:,0,0], 'k', label='daily mean')
        axs[1].plot(smooth[:,0,0], 'b', label='smoothed')
        axs[1].legend()
        axs[1].set_ylim(axs[0].get_ylim())
        axs[1].set_xlabel('time [d]')
        axs[1].set_ylabel(name)
        plt.show()

    return smooth

def raw_precip_to_spinup(precip_rain, precip_snow, air_temp_spinup_K, plot=False):
    """Stacks years and generates a total averaged precip rate for snow and rain."""
    precip = smooth_raw_to_spinup(precip_rain + precip_snow, 'total precip', plot)

    assert(air_temp_spinup_K.shape == precip.shape)
    spinup_rain = np.where(air_temp_spinup_K >= 273.15, precip, 0.)
    spinup_snow = np.where(air_temp_spinup_K < 273.15, precip, 0.)
    return spinup_rain, spinup_snow

def raw_to_spinup(raw, spinup_n_years, plot=False):
    """From raw years, generate a single year of spinup data."""
    spinup = dict()
    all_keys = list(raw.keys())

    # copy over x and y from raw, these won't change
    for key in _prewritten:
        if key.startswith('time'):
            # time is simply daily data
            if '[s]' in key:
                spinup[key] = 86400. * np.arange(0., spinup_n_years * 365., 1.)
            elif '[d]' in key:
                spinup[key] = np.arange(0., spinup_n_years * 365., 1.)
        else:
            spinup[key] = raw[key]
        all_keys.remove(key)

    # smooth all other keys
    for key in all_keys:
        if not key.startswith('precipitation'):
            spinup[key] = smooth_raw_to_spinup(raw[key], key, plot)

    # precip requires partitioning
    rain, snow = raw_precip_to_spinup(raw['precipitation rain [m s^-1]'],
                                      raw['precipitation snow [m SWE s^-1]'],
                                      spinup['air temperature [K]'][:], plot)
    spinup['precipitation rain [m s^-1]'] = rain
    spinup['precipitation snow [m SWE s^-1]'] = snow

    if plot:
        fig, ax = plt.subplots(1,1)
        ax.plot(rain[:,0,0] * 86400 * 1000, 'b', label='rain')
        ax.plot(snow[:,0,0] * 86400 * 1000, 'c', label='snow')
        ax.set_xlabel('time [d]')
        ax.set_ylabel('rate [mm / day]')
        ax.legend()
        plt.show()

    return spinup

def write_ats(dat, attrs, filename):
    """Accepts a dictionary of ATS data and writes it to HDF5 file."""
    logging.info('Writing ATS file: {}'.format(filename))
    ntimes = len(dat['time [s]'])
    
    with h5py.File(filename, 'w') as fid:
        for key in _prewritten:
            fid.create_dataset(key, data=dat[key])

        for key, val in dat.items():
            if not key in _prewritten:
                group = fid.create_group(key)
                for i in range(ntimes):
                    doy = i%365
                    group.create_dataset(str(i), data=dat[key][doy,:,:])

        for key, val in attrs.items():
            fid.attrs[key] = val
    return


def read_raw(fid):
    raw = dict()
    for key in _prewritten:
        raw[key] = fid[key][:]

    ntimes = len(raw['time [s]'])
    for key in fid.keys():
        if not key in raw:
            d = np.array([fid[key][str(i)][:] for i in range(ntimes)])
            raw[key] = d
            assert(d.shape[0] == ntimes)
    return raw
                    
                
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('input_filename', help='Input filename, an ATS-style HDF5 file.')
    # parser.add_argument('-s', '--spinup', action='store_true',
    #                     help='Write a spinup file, which is ultra smoothed.')
    parser.add_argument('-n', '--num-years', type=int, default=1,
                        help='Number of years to write.')

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-d', '--directory', default='.',
                        help='Directory in which to place output files, with default names.')
    group.add_argument('-o', '--output-filename',
                        help='Full path to output filename')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='Plot data.')


    args = parser.parse_args()
    if args.output_filename is None:
        if args.input_filename.startswith('daymet_raw'):
            args.output_filename = os.path.join(args.directory,
                                                args.input_filename.replace('daymet_raw', 'spinup', 1))
        else:
            args.output_filename = os.path.join(args.directory,
                                                'spinup_'+args.input_filename)

    # read the data, create the spinup data
    with h5py.File(args.input_filename, 'r') as fid:
        raw = read_raw(fid)
        attrs = dict(fid.attrs)

    spinup = raw_to_spinup(raw, args.num_years, args.plot)

    # write to disk
    write_ats(spinup, attrs, args.output_filename)

        
        
    
    
