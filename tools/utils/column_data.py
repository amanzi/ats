#!/usr/bin/env python

"""Loads, sorts, writes, and plots column data from a 1D (potentially disordered) simulation into a H5 file.

Usage: column_data.py

Loads ATS simulation data in the current directory and writes a sorted
file for pressure and temperature named column_data.h5 which can then
be read by other simulations looking to interpolate a 1D initial
condition.

NOTE: this is deprecated, please prefer to use plot_column_data.py.  However,
it was used in a lot of spinup setups, so is kept for posterity.

"""

import sys,os
import numpy as np
import h5py
import mesh

def meshX(coord=2, filename="visdump_mesh.h5", directory=".", key=None):
    return mesh.meshElemCentroids(filename,directory,key)[:,coord]

def fullname(varname):
    fullname = varname
    if not '.cell.' in fullname:
        fullname = fullname+'.cell.0'
    return fullname

def column_data(varnames, keys='all', directory=".", filename="visdump_data.h5",
                mesh_filename="visdump_mesh.h5", coord=2, deformable=False):
    """Returns data of shape ( len(varnames+1), len(keys), n_cells )"""
    if type(varnames) is str:
        varnames = [varnames,]

    z = meshX(coord, mesh_filename, directory)

    with h5py.File(os.path.join(directory,filename),'r') as dat:
        keys_avail = dat[fullname(varnames[0])].keys()
        keys_avail.sort(lambda a,b: int.__cmp__(int(a),int(b)))

        if keys == 'all':
            keys = keys_avail
        elif type(keys) is str:
            keys = [keys,]
        elif type(keys) is int:
            keys = [keys_avail[keys],]
        elif type(keys) is slice:
            keys = keys_avail[keys]

        vals = np.zeros((len(varnames)+1, len(keys), len(z)), 'd')
        for i,key in enumerate(keys):
            if deformable:
                z = meshX(coord, mesh_filename, directory, key)
            vals[0,i,:] = z
            for j,varname in enumerate(varnames):
                vals[j+1,i,:] = dat[fullname(varname)][key][:,0]

    # sort in z coordinate
    return vals[:,:,vals[0,0,:].argsort()]

def getFigs(inset, is_temp, figsize=(12,3)):
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=figsize)
    axs = []
    if is_temp:
        axs.append(fig.add_subplot(131))
        axs.append(fig.add_subplot(132))
        axs.append(fig.add_subplot(133))
        if inset:
            axs.append(fig.add_axes([0.6,0.5,0.25, 0.25]))
    else:
        axs.append(fig.add_subplot(121))
        axs.append(fig.add_subplot(122))
        if inset:
            axs.append(fig.add_axes([0.86,0.67,0.1, 0.25]))
    return fig,axs


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate columnar data from an unstructured column run.')
    parser.add_argument('-t', '--temperature', action='store_true',
                        help='include temperature data')
    parser.add_argument('z0', metavar='TOP_SURFACE_ELEVATION', type=float,
                        help='elevation of the top surface of the column run')
    options = parser.parse_args()
    
    with h5py.File("column_data.h5", 'w') as fout:
        to_read = ['pressure']
        if options.temperature:
            to_read.append("temperature")
        dat = column_data(to_read, keys=-1)
        z_depth = options.z0 - dat[0,0,:] # correction to get surface set, depth coordinate
        fout.create_dataset("z", data=np.flipud(z_depth))
        fout.create_dataset("pressure", data=np.flipud(dat[1,0,:]))
        if options.temperature:
            fout.create_dataset("temperature", data=np.flipud(dat[2,0,:]))
