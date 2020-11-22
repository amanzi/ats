#!/usr/bin/env python
"""Loads and plots timestep history for a given run."""
from __future__ import print_function
from __future__ import division

import numpy as np

def print_headers():
    print("cycle, time, dt, iteration count, wallclock avg (s)")

def parse_logfile(fid, wallclock=False):
    """Reads a file, and returns a list of good and bad timesteps.

    Each are a list of 3-tuples: (step number, time of step, dt)
    """
    
    data = []
    faildata = []

    # find number of "cycles"
    total_cycles = 0
    for line in fid:
        if "Cycle =" in line:
            sline = line.split()
            cyc = int(sline[4][:-1])
            time = float(sline[8][:-1])
            dt = float(sline[12])

            if len(data) > 0 and data[-1][0] == cyc:
                faildata.append(data.pop())
            data.append([cyc,time,dt])

    return np.array(data), np.array(faildata)


def get_axs():
    """Gets a figure and list of axes for plotting."""
    return plt.subplots(1,2)

def decorate_axs(axs):
    """Adds legends, labels, limits."""
    axs[0].set_xlabel("time [days]")
    axs[0].set_ylabel("dt [days]")
    axs[0].legend(loc='lower left')
    axs[1].set_xlabel("cycles [-]")
    axs[1].set_ylabel("dt [days]")
    return

def plot(data, axs, color, label, symbol='x'):
    """Plot the data."""
    axs[0].semilogy(data[0][:,1], data[0][:,2], '-'+symbol, color=color, label=label)
    if data[1].shape != (0,):
        axs[0].semilogy(data[1][:,1], data[1][:,2], symbol, color=color)
    axs[1].semilogy(data[0][:,0], data[0][:,2], '-'+symbol, color=color)
    if data[1].shape != (0,):
        axs[1].semilogy(data[1][:,0], data[1][:,2], symbol, color=color)

def write_to_file(data, fnamebase):
    """Writes the data to a file for future reading"""
    fname = fnamebase+".npz"
    np.savez(fname, good_timesteps=data[0], bad_timesteps=data[1])

def read_from_file(fname):
    """Reads a .npz file"""
    read = np.load(fname)
    return [read["good_timesteps"], read["bad_timesteps"]]
            
if __name__ == "__main__":
    import sys,os
    from matplotlib import pyplot as plt

    import argparse
    import colors
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("LOG_FILES", nargs="+", type=str,
                        help="List of logfiles to parse.")
    parser.add_argument("--colors", "-c", type=colors.float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    parser.add_argument("--overwrite", "-o", action="store_true",
                        help="Do not use any existing .npz file -- instead reload from the logfile and overwrite the .npz file.")
    args = parser.parse_args()

    if args.colors is not None:
        cm = colors.cm_mapper(0,1,args.colormap)
        args.colors = [cm(c) for c in args.colors]
    else:
        args.colors = colors.enumerated_colors(len(args.LOG_FILES))
        
    fig, axs = get_axs()
        
    for fname, color in zip(args.LOG_FILES, args.colors):
        if fname.endswith(".npz"):
            data = read_from_file(fname)
        elif os.path.isfile(fname+".npz") and not args.overwrite:
            data = read_from_file(fname+".npz")
        else:
            with open(fname,'r') as fid:
                data = parse_logfile(fid)
            write_to_file(data, fname)
                
        plot(data, axs, color, fname)

    decorate_axs(axs)
    plt.show()
    sys.exit(0)


