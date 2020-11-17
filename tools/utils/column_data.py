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
import ats_xdmf

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate columnar data from an unstructured column run.')
    parser.add_argument('-t', '--temperature', action='store_true',
                        help='include temperature data')
    parser.add_argument('z0', metavar='TOP_SURFACE_ELEVATION', type=float,
                        help='elevation of the top surface of the column run')
    options = parser.parse_args()

    vis = ats_xdmf.VisFile()
    vis.loadMesh(columnar=True)
    vis.filterIndices(-1)

    pres = vis.getArray('pressure')
    if options.temperature:
        temp = vis.getArray('temperature')

    zc = vis.centroids[:,2]
    z_depth = options.z0 -  zc

    with h5py.File("column_data.h5", 'w') as fout:
        fout.create_dataset('z', data=np.flipud(z_depth))
        fout.create_dataset('pressure', data=np.flipud(pres[0,:]))
        if options.temperature:
            fout.create_dataset('temperature', data=np.flipud(temp[0,:]))

