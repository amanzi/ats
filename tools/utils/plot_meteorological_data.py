#!/usr/bin/env python
"""Plots meteorological data, using expected default names, from an ATS h5 input file.

Usage: plot_meteorological_data.py INPUT.h5
"""

import h5py
import numpy as np
from matplotlib import pyplot as plt

import colors


names1 = ["air temperature [K]",
         "relative humidity [-]",
         "incoming shortwave radiation [W m^-2]",
         "incoming longwave radiation [W m^-2]",
         "wind speed [m s^-1]"]

precip1 = ["precipitation rain [m s^-1]",
          "precipitation snow [m SWE s^-1]"]


names2 = ["Ta", "RH", "Qswin", "Qlwin", "Us"]
precip2 = ["Pr", "Ps"]

def plot_met(fname, axs, color='b', end_time_in_years=np.inf, style='-'):
    axs = axs.ravel()

    with h5py.File(fname, 'r') as fid:
        try:
            time = fid['time [s]'][:]/86400.0/365.0
        except KeyError:
            try:
                time = fid['time'][:]/86400.0/365.0
            except KeyError:
                raise KeyError('Missing time entry "time [s]"')
            else:
                names = names2
                precip = precip2
        else:
            names = names1
            precip = precip1

        time = np.squeeze(time)            

        if end_time_in_years > time[-1]:
            end = len(time)
        else:
            end = np.where(time > end_time_in_years)[0][0]
        
        for i,n in enumerate(names):
            if n in fid.keys():
                axs[i].plot(time[0:end], np.squeeze(fid[n][0:end]), style, color=color)
            axs[i].set_title(n)

        for p, s2 in zip(precip, ['-','--']):
            if p in fid.keys():
                axs[i+1].plot(time[0:end], np.squeeze(fid[p][0:end]), s2, color=color)
        axs[i+1].set_title("precip (solid=rain, dash=snow) [m SWE s^-1]")


def get_axs():
    return plt.subplots(3,2, figsize=(10,10))
        
if __name__ == "__main__":
    import sys
    import colors
    import argparse

    parser = argparse.ArgumentParser(description="Plot met data from an ATS input h5 file using default names.")
    parser.add_argument("INPUT_FILES", nargs="+", type=str,
                        help="List of logfiles to parse.")
    args = parser.parse_args()

    color_list = colors.enumerated_colors(len(args.INPUT_FILES))

    fig, axs = get_axs()

    for fname, color in zip(args.INPUT_FILES, color_list):
        plot_met(fname, axs, color)

    plt.show()
    sys.exit(0)
