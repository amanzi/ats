"""Functions for parsing and working with ATS meshes."""

import sys,os
import numpy as np
import h5py

import geometry

elem_type = {5:'QUAD',
             8:'PRISM',
             9:'HEX',
             4:'TRIANGLE'
             }

def meshElemXYZ(filename="ats_vis_mesh.h5", directory=".", key="/0"):
    """Returns a set of coordinates of nodes in the mesh, along with a elemental conection list."""

    with h5py.File(os.path.join(directory, filename), 'r') as dat:
        #if key is None:
            #key = dat.keys()[0]
        
        key="/0"
        mesh = dat[key]['Mesh']
        elem_conn = mesh['MixedElements'][:,0]

        etype = elem_type[elem_conn[0]]
        if (etype == 'PRISM'):
            nnodes_per_elem = 6
            nnodes_per_face = 3
        elif (etype == 'HEX'):
            nnodes_per_elem = 8
            nnodes_per_face = 4
        elif (etype == 'QUAD'):
            nnodes_per_elem = 4
            nnodes_per_face = 2
        elif (etype == 'TRIANGLE'):
            nnodes_per_elem = 3
            nnodes_per_face = 2

        n_elems = len(elem_conn) / (nnodes_per_elem+1)
        coords = dict(zip(mesh['NodeMap'][:,0], mesh['Nodes'][:]))

    conn = elem_conn.reshape((int(n_elems), nnodes_per_elem+1))
    assert np.all(conn[:,0] == elem_conn[0])
    return etype, coords, conn


def meshElemCentroids(filename="ats_vis_mesh.h5", directory=".", key=None):
    """Gets centroids of cells in the mesh."""
    
    etype, coords, conn = meshElemXYZ(filename,directory,key)

    centroids = np.zeros((len(conn),3),'d')
    for i,elem in enumerate(conn):
        elem_coords = np.array([coords[gid] for gid in elem[1:]])
        elem_z = np.mean(elem_coords, axis=0)
        centroids[i,:] = elem_z
    return centroids
    
def meshCentroidsStructuredOrdering_3D(order=["x","z"], filename="ats_vis_mesh.h5", directory="."):
    """Returns centroids of a structured mesh ordered in the natural ordering.
    order: a list given to numpy sort(), headings are "x","y","z".
           ["x","z"] (default) gives coordinates in z-fastest
    """
    
    c = meshElemCentroids(filename, directory)

    coords_a = np.array([(i,c[i,0],c[i,1],c[i,2]) for i in range(c.shape[0])],
                        dtype=[('id',int),('x',float),('y',float),('z',float)])
    coords_a.sort(order=order)
    return coords_a
    

def meshElemVolume(filename="ats_vis_mesh.h5", directory="."):
    """Returns the volume of each cell.
    Note this is the cell volume in 3D or the face area in 2D surface mesh."""
    
    etype, coords, conn = meshElemXYZ(filename, directory)
    volumes = np.zeros((len(conn),),'d')

    if etype == 'PRISM' or etype == 'HEX':
        raise NotImplementedError("meshElemVolume not implemented for 3D.")

    elif etype == 'QUAD' or etype == 'TRIANGLE':
        for i,elem in enumerate(conn):
            elem_coords = [coords[gid] for gid in elem[1:]]
            volumes[i] = geometry.poly_area(elem_coords)
            
    else:
        raise NotImplementedError("meshElemVolume not implemented for element type: %s"%etype)
        
    return volumes