"""Loads and/or plots 2D, topologlically structured data on quadrilaterals using matplotlib.
"""

import sys,os
import numpy as np
import h5py
import mesh
import colors

def fullname(varname):
    fullname = varname
    if not '.cell.' in fullname:
        fullname = fullname+'.cell.0'
    return fullname


def transect_data(varnames, keys='all', directory=".", filename="visdump_data.h5",
                  mesh_filename="visdump_mesh.h5", coord_order=None, deformable=False, return_map=False):
    """Pulls simulation output into structured 2D arrays for transect-based, (i,j) indexing.

    Input:
      varnames       | A list of variable names to pull, e.g.
                     |  ['saturation_liquid', 'saturation_ice'], or a single variable
                     |  name, e.g. 'saturation_liquid'
      keys           | Indices of timesteps to pull.  Either an int (i.e. 0, -1, etc) 
                     |  for the kth timestep, or a list of ints, or 'all'.
      directory      | Directory of the run.  Defaults to '.'
      filename       | Filename of the run.  Defaults to 'visdump_data.h5'
      mesh_filename  | Filename of the mesh.  Defaults to 'visdump_mesh.h5'
      coord_order    | Order of the transect coordinates.  Defaults to ['x','z'].  The 
                     |  mesh is sorted in this order.
      deformable     | Is the mesh deforming?
      return_map     | See return value below.
   
    Output:
      Output is an array of shape:
      ( len(varnames+2), len(keys), n_cells_coord_order[0], n_cells_coord_order[1] )
    
      data[0,0,:,:] is the coord_order[0] centroid
      data[1,0,:,:] is the coord_order[1] centroid
      data[i+2,k,:,:] is the ith varname data at the kth requested timestep, sorted in 
                      the same way as the centroids.

      Note that the data is re-ordered in INCREASING coordinate, i.e. bottom to top in z.

      If return_map is True, then returns a tuple, (data, map) where
      map is a (NX,NZ) array of integers specifying which global id
      corresponds to the (i,j) cell.  This is useful for mapping input
      data back INTO the unstructured mesh.

    Example usage:  
      Calculate and plot the thaw depth at step 5.

      // Pull saturation ice -- TD is where sat ice = 0."
      data = transect_data(['saturation_ice', 5)

      // x coordinate for plotting
      x = data[0,0,:,0]
    
      // for each column, find highest z where sat_ice > 0.
      td_i = np.array([np.where(data[2,0,i,:] > 0.)[0][-1] for i in range(data.shape[2])])

      // now that we have an index into the highest cell with ice, determine td as the 
      // mean of the highest cell with ice and the one above that.  Note this assumes
      // all columns have some thawing.
      td_z = np.array( [  (dat[1,0,i,td_i[i]] + dat[1,0,i,td_i[i+1]]) / 2. 
                              for i in range(len(td_i)) ] )

      plt.plot(x, td_z)

    """
    if coord_order is None:
        coord_order = ['x','z']

    if type(varnames) is str:
        varnames = [varnames,]

    # get centroids
    xyz = mesh.meshElemCentroids(mesh_filename, directory)

    # round to avoid issues
    xyz = np.round(xyz, decimals=5)

    # get ordering of centroids
    dtype = [(coord_order[0], float), (coord_order[1], float)]
    num_order = []
    for i in coord_order:
        if i == 'x':
            num_order.append(0)
        elif i == 'y':
            num_order.append(1)
        elif i == 'z':
            num_order.append(2)

    xyz_sort_order = np.array([tuple([xyz[i,x] for x in num_order]) for i in range(len(xyz))], dtype=dtype)
    xyz_sorting = xyz_sort_order.argsort(order=coord_order)

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
        elif type(keys) is list:
            if all(type(k) is int for k in keys):
                keys = [keys_avail[k] for k in keys]
            elif all(type(k) is str for k in keys):
                pass
            else:
                raise RuntimeError("Keys requested cannot be processed -- should be 'all', int, or str key, or list of ints or strs.")
                

        # get data
        vals = np.zeros((len(varnames)+2, len(keys), len(xyz)), 'd')

        for i,key in enumerate(keys):
            if deformable:
                xyz = mesh.meshElemCentroids(mesh_filename, directory)
            vals[0,i,:] = xyz[xyz_sorting,num_order[0]]
            vals[1,i,:] = xyz[xyz_sorting,num_order[1]]
            for j,varname in enumerate(varnames):
                vals[j+2,i,:] = dat[fullname(varname)][key][:,0][xyz_sorting]

    # reshape the data
    # determine nx
    nx = len(set(vals[0,0,:]))
    nz = vals.shape[2] / nx
    if (nx * nz != vals.shape[2]):
        raise RuntimeError("Assumption about first coordinate being cleanly binnable is falling apart -- ask Ethan to rethink this algorithm!")
    shp = vals.shape

    if not return_map:
        return vals.reshape(shp[0], shp[1], nx, nz)
    else:
        return vals.reshape(shp[0], shp[1], nx, nz), xyz_sorting.reshape(nx, nz)


def plot(dataset, ax, cax=None, vmin=None, vmax=None, cmap="jet",
         label=None, mesh_filename="visdump_mesh.h5", directory=".", y_coord=0.0,
         linewidths=1):
    """Draws a dataset on an ax."""
    import matplotlib.collections
    from matplotlib import pyplot as plt

    if vmin is None:
        vmin = dataset.min()
    if vmax is None:
        vmax = dataset.max()

    # get the mesh and collapse to 2D
    etype, coords, conn = mesh.meshElemXYZ(filename=mesh_filename, directory=directory)
    if etype is not 'HEX':
        raise RuntimeError("Only works for Hexs")

    coords2 = np.array([[coords[i][0::2] for i in c[1:] if abs(coords[i][1] - y_coord) < 1.e-8] for c in conn])
    try:
        assert coords2.shape[2] == 2
        assert coords2.shape[1] == 4
    except AssertionError:
        print(coords2.shape)
        for c in conn:
            if len(c) != 9:
                print c
                raise RuntimeError("what is a conn?")
            coords3 = np.array([coords[i][:] for i in c[1:] if abs(coords[i][1] - y_coord) < 1.e-8])
            if coords3.shape[0] != 4:
                print coords
                raise RuntimeError("Unable to squash to 2D")

    # reorder anti-clockwise
    for i,c in enumerate(coords2):
        centroid = c.mean(axis=0)
        def angle(p1,p2):
            a1 = np.arctan2((p1[1]-centroid[1]),(p1[0]-centroid[0]))
            a2 = np.arctan2((p2[1]-centroid[1]),(p2[0]-centroid[0]))
            if a1 < a2:
                return -1
            elif a2 < a1:
                return 1
            else:
                return 0

        c2 = np.array(sorted(c,angle))
        coords2[i] = c2

    polygons = matplotlib.collections.PolyCollection(coords2, edgecolor='k', cmap=cmap, linewidths=linewidths)
    polygons.set_array(dataset)
    polygons.set_clim(vmin,vmax)
    ax.add_collection(polygons)

    xmin = min(c[0] for c in coords.itervalues())
    xmax = max(c[0] for c in coords.itervalues())
    zmin = min(c[2] for c in coords.itervalues())
    zmax = max(c[2] for c in coords.itervalues())

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(zmin,zmax)

    if cax is not None:
        cb = plt.colorbar(polygons, cax=cax)
        if label is not None:
            cb.set_label(label)
        
    return ((xmin,xmax),(zmin,zmax))
    
