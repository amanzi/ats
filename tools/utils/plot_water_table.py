#!/usr/bin/env python
"""
Plot water table as a function of time.

This is currently only useful on a single, 1D column.
"""

import h5py
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm
import ats_xdmf, column_data
    

def water_table_no_surface(dirname, datum=None, patm=101325.):
    # get the columnar pressure
    col_dat = column_data.column_data(['pressure',], directory=dirname)
    wt = np.zeros((col_dat.shape[1],),'d')
    if datum is None:
        datum = column_data[0,0,-1]
    
    for k in range(col_dat.shape[1]):
        if col_dat[1,k,-1] > patm:
            wt[k] = (col_dat[1,k,-1]-patm) / 1000. / 9.81 - datum  # hydrostatic above
        else:
            wt_index = np.where(col_dat[1,k,:] >= patm)[0]
            if len(wt_index) is 0:
                wt_index = 0
            else:
                wt_index = wt_index[-1]
            wt[k] = col_dat[0,k,wt_index] - datum
    return times, wt



def water_table(dirname, datum=None, v86=False, patm=101325.):
    prefix = "" if v86 else "surface-"

    # get the ponded depth
    vf = ats_xdmf.VisFile(dirname, filename="visdump_surface_data.h5")
    keys = vf.cycles
    times = vf.times
    pd_key = prefix + "ponded_depth.cell.0"
    pd = np.array([vf.d[pd_key][k][0] for k in keys])
    elev_surf = vf.d[prefix + "elevation.cell.0"][keys[0]][0]
    vf.close()
    
    if datum is None:
        datum = elev_surf
    datum_offset = elev_surf - datum

    # get the columnar pressure
    wt = np.zeros((len(keys),),'d')
    col_dat = column_data.column_data(['pressure',], directory=dirname)
    for k in range(col_dat.shape[1]):
        if pd[k] > 0:
            wt[k] = pd[k] + datum_offset
        else:
            wt_index = np.where(col_dat[1,k,:] >= patm)[0]
            if len(wt_index) is 0:
                wt_index = 0
            else:
                wt_index = wt_index[-1]
            wt[k] = col_dat[0,k,wt_index] - datum
    return times, wt


def plot_water_table(times, wt, ax, **kwargs):
    ax.plot(times, wt, **kwargs)

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
    parser.add_argument("-p", "--permafrost", action="store_true",
                        help="Permafrost run, plot ice table too.")
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

        times,wt = water_table(fname)
        plot_water_table(times, wt, axs, color=c, label=fname)

        if args.permafrost:
            import plot_thaw_depth
            times, td = plot_thaw_depth.thaw_depth(fname)
        plot_thaw_depth.plot_thaw_depth(times, td, axs, linestyle='--', color=c)

    axs.plot([times[0], times[-1]], [0.,0.], '--', color='k')
    plt.tight_layout()
    axs.legend()
    plt.show()
