#!/usr/bin/env python
"""Loads and plots timestep history for a given run."""
import numpy as np
from matplotlib import pyplot as plt

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

    ngood = len(data)
    nbad = len(faildata)
    print(f"  total timesteps = {ngood+nbad}")
    print(f"  bad timesteps = {nbad}")
    return np.array(data), np.array(faildata)


def get_axs():
    """Gets a figure and list of axes for plotting."""
    return plt.subplots(1,2)

def decorate_axs(axs, units='d'):
    """Adds legends, labels, limits."""
    
    axs[0].set_xlabel(f"time [{units}]")
    axs[0].set_ylabel("dt [days]")
    axs[0].legend(loc='lower left')
    axs[1].set_xlabel("cycles [-]")
    axs[1].set_ylabel(f"dt [days]")
    return

def plot(data, axs, color, label, symbol='x', units='s'):
    """Plot the data."""
    if units == 's':
        factor = 86400
    elif units == 'd':
        factor = 1
    elif units == 'y':
        factor = 1.0 / 365
    
    axs[0].semilogy(data[0][:,1]*factor, data[0][:,2], '-'+symbol, color=color, label=label)
    if data[1].shape != (0,):
        axs[0].semilogy(data[1][:,1]*factor, data[1][:,2], symbol, color=color)

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
    ngood = len(read['good_timesteps'])
    nbad = len(read['bad_timesteps'])
    print(f"File: {filename}")
    print(f"  total timesteps = {ngood+nbad}")
    print(f"  bad timesteps = {nbad}")
    return [read["good_timesteps"], read["bad_timesteps"]]
            
if __name__ == "__main__":
    import sys,os
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
    parser.add_argument("--units", type=str, default='d', choices=['y', 'd', 's'],
                        help="Set the units of the x axis -- one of 'y', 'd', 's'")
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
            print(f"File: {fname}")
            with open(fname,'r') as fid:
                data = parse_logfile(fid)
            write_to_file(data, fname)
                
        plot(data, axs, color, fname, units=args.units)

    decorate_axs(axs, args.units)
    plt.show()
    sys.exit(0)


