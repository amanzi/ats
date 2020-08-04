#!/bin/env python

import numpy as np
from matplotlib import pyplot as plt
import argparse
import sys
import parse_ats

def load(fname, density):
    dat = np.loadtxt(fname) # units s, mol/s
    dat[:,0] = dat[:,0] / 86400. # convert to days
    dat[:,1] = dat[:,1] / density * 86400 # convert to m^3/d
    return dat

def plot(data, format='-', color='b', name=None, ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.plot(data[:,0], data[:,1], format, color=color, label=name)
    ax.set_xlabel("time [days]")
    ax.set_ylabel("runoff [m^3 / day]")
    return ax

def load_area_rain(args):
    k,t,d = parse_ats.readATS(args.directory, args.filename)
    cv = d[args.area_key][k[0]][:]
    area = cv.sum()
    rain = np.array([(d[args.rainfall_rate_key][key][:] * cv).sum() for key in k]) * 86400
    return area, t*365.25, rain # units m^2, days, m^3/s

def plot_rain(area, t, rain, format='--', color='k', ax=None):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    ax.plot(t, rain, format, color=color, label="rainfall rate")
    return ax

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Plot discharge observation from ATS run")
    parser.add_argument("runoff_filename", type=str, help="Runoff observation filename.")
    parser.add_argument("-p", "--plot-rainfall", action="store_true", help="Plot rainfall rate as an asymptotic limit.")
    parser.add_argument("-d", "--directory", type=str, help="Simulation output directory", default='.')
    parser.add_argument("-f", "--filename", type=str, help="Simulation surface output filename", default="visdump_surface_data.h5")
    parser.add_argument("-r", "--rainfall-rate-key", type=str, help="Rainfall rate variable name", default="surface-mass_source.cell.0")
    parser.add_argument("-a", "--area-key", type=str, help="Surface cell area variable name", default="surface-cell_volume.cell.0")
    parser.add_argument("--density", type=float, help="Density of water", default=55000.)
    args = parser.parse_args()

    ax = None
    if args.plot_rainfall:
        area, time, rain = load_area_rain(args)
        ax = plot_rain(area, time, rain)

    plot(load(args.runoff_filename, args.density), ax=ax)
    plt.show()
    sys.exit(0)
