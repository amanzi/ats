"""Functionality to use daily meterological data to form smoother daily data."""

import sys, os
import numpy as np
import scipy.signal
import h5py
from matplotlib import pyplot as plt
import logging

def stack_years(raw):
    """Checks that data can be cast in years and stacks it into a 365 x NYEAR array."""
    if len(raw)%365 != 0:
        raise RuntimeError('Raw data provided is not module 365 (len = %d)'%len(raw))

    nyears = int(np.round(len(raw) / 365))
    raw = raw.reshape((nyears,365)).transpose()
    return nyears, raw

def smooth_raw_to_spinup(raw, name=None, plot=False):
    """Stacks years and generates a single, smooth year.

    Assumes input met data is daily data, and in 365-day year format such that 
    len(met)%365 == 0.
    """
    nyears, raw = stack_years(raw)
    
    mean = raw.mean(axis=1)
    smooth = scipy.signal.savgol_filter(mean, 61, 2, mode='wrap')

    if plot:
        fig, axs = plt.subplots(2,1)
        for i in range(nyears):
            axs[0].plot(raw[:,i])
        axs[0].set_xticklabels(list())
        axs[1].plot(mean, 'k', label='daily mean')
        axs[1].plot(smooth, 'b', label='smoothed')
        axs[1].legend()
        axs[1].set_ylim(axs[0].get_ylim())
        axs[1].set_xlabel('time [d]')
        axs[1].set_ylabel(name)
        plt.show()

    return smooth

def raw_precip_to_spinup(precip_rain, precip_snow, air_temp_spinup_K, plot=False):
    """Stacks years and generates a single averaged precip rate for snow and rain."""
    nyears, precip = stack_years(precip_rain+precip_snow)
    mean_rate = precip.mean()

    spinup_rain = np.where(air_temp_spinup_K >= 273.15, mean_rate, 0.)
    spinup_snow = np.where(air_temp_spinup_K < 273.15, mean_rate, 0.)

    if plot:
        fig, ax = plt.subplots(1,1)
        ax.plot(spinup_rain * 86400 * 1000, 'b', label='rain')
        ax.plot(spinup_snow * 86400 * 1000, 'c', label='snow')
        ax.set_xlabel('time [d]')
        ax.set_ylabel('rate [mm / day]')
        ax.legend()
        plt.show()

    return spinup_rain, spinup_snow

def raw_to_spinup(raw, spinup_n_years, plot=False):
    """From raw years, generate a single year of spinup data."""
    spinup = dict()
    
    # air temp, relative humidity, wind speed, and incoming shortwave are simply smoothed
    for k in ['air temperature [K]', 'wind speed [m s^-1]',
              'relative humidity [-]', 'incoming shortwave radiation [W m^-2]']:
        spinup[k] = smooth_raw_to_spinup(raw[k][:], k, plot)

    # precip requires partitioning
    rain, snow = raw_precip_to_spinup(raw['precipitation rain [m s^-1]'][:],
                                      raw['precipitation snow [m SWE s^-1]'][:],
                                      spinup['air temperature [K]'][:], plot)
    spinup['precipitation rain [m s^-1]'] = rain
    spinup['precipitation snow [m SWE s^-1]'] = snow

    # tile all data to repeat n_year times
    for key, val in spinup.items():
        spinup[key] = np.tile(spinup[key], spinup_n_years)
    
    # time is simply daily data
    spinup['time [s]'] = 86400. * np.arange(0., spinup_n_years * 365., 1.)
    return spinup

def write_ats(dat, attrs, filename):
    """Accepts a dictionary of ATS data and writes it to HDF5 file."""
    logging.info('Writing ATS file: {}'.format(filename))
    with h5py.File(filename, 'w') as fid:
        for key, val in dat.items():
            fid.create_dataset(key, data=val)

        for key, val in attrs.items():
            fid.attrs[key] = val
    return


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
        spinup = raw_to_spinup(fid, args.num_years, args.plot)
        attrs = dict(fid.attrs)

    # write to disk
    write_ats(spinup, attrs, args.output_filename)

        
        
    
    
