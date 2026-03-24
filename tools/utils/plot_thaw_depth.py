#!/usr/bin/env python
"""
Plot thaw depth as a function of time.

This is currently only useful on a single, 1D column.
"""

import h5py
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm
import ats_xdmf, column_data
    

def thaw_depth(dirname, datum=None, v86=False, T0=273.15):
    prefix = "" if v86 else "surface-"

    # get the unfrozen fraction, ponded depth
    vf = ats_xdmf.VisFile(dirname, filename="visdump_surface_data.h5")
    keys = vf.cycles
    times = vf.times
    pd_key = prefix + "ponded_depth.cell.0"
    uf_key = prefix + "unfrozen_fraction.cell.0"
    pd = np.array([vf.d[pd_key][k][0] for k in keys])
    uf = np.array([vf.d[uf_key][k][0] for k in keys])
    elev_surf = vf.d[prefix + "elevation.cell.0"][keys[0]][0]
    vf.close()
    
    if datum is None:
        datum = elev_surf
    datum_offset = elev_surf - datum

    td = np.zeros((len(keys),),'d')
    # get the columnar temperature
    col_dat = column_data.column_data(['temperature',], directory=dirname)
    for k in range(col_dat.shape[1]):
        if uf[k] == 0.:
            td[k] = pd[k] + datum_offset
        else:
            td[k] = pd[k] * (1 - uf[k]) + datum_offset

            still_unfrozen = True
            c = col_dat.shape[2]-1
            while still_unfrozen and c >= 0:
                if col_dat[1,k,c] > T0:
                    td[k] = col_dat[0,k,c] - datum
                    c = c - 1
                else:
                    still_unfrozen = False
    return times, td

def plot_thaw_depth(times, td, ax, **kwargs):
    ax.plot(times, td, **kwargs)

def get_axs():
    f,axs = plt.subplots(1,1)
    return f,axs
        
if __name__ == "__main__":
    import sys

    import argparse
    import colors
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("INPUT_DIRECTORIES", nargs="+", type=str,
                        help="List of directories to plot.")
    parser.add_argument("--colors", "-c", type=colors.float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    args = parser.parse_args()

    cm = colors.cm_mapper(0,1,args.colormap)
    fig, axs = get_axs()

    fnames = args.INPUT_DIRECTORIES
    for i,fname in enumerate(fnames):
        if args.colors is None:
            if len(fnames) > 1:
                c = cm(float(i)/(len(fnames)-1))
            else:
                c = 'b'
        else:
            if type(args.colors[i]) is float:
                c = cm(args.colors[i])
            else:
                c = args.colors[i]

        times,td = thaw_depth(fname)
        plot_thaw_depth(times, td, axs, color=c, label=fname)

    plt.tight_layout()
    axs.legend()
    plt.show()
