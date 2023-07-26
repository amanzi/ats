
"""Functions for parsing Amanzi/ATS XDMF visualization files."""
import sys,os
import numpy as np
import h5py
import matplotlib.collections

def valid_data_filename(domain, format=None):
    """The filename for an HDF5 data filename formatter"""
    if format is None:
        format = 'ats_vis_{}_data.h5'

    if domain == 'domain':
        domain = ''
    fname = format.format(domain)
    fname = fname.replace('__', '_')
    return fname
    
def valid_mesh_filename(domain, format=None):
    """Argparse validator for an HDF5 mesh filename formatter"""
    if format is None:
        format = 'ats_vis_{}_mesh.h5'
    return valid_data_filename(domain, format)


class VisFile:
    """Class managing the reading of ATS visualization files."""
    def __init__(self, directory='.', domain=None, filename=None, mesh_filename=None, time_unit='yr'):
        """Create a VisFile object.

        Parameters
        ----------
        directory : str, optional
          Directory containing vis files.  Default is '.'
        domain : str, optional
          Amanzi/ATS domain name.  Useful in variable names, filenames, and more.
        filename : str, optional
          Filename of h5 vis file.  Default is 'ats_vis_DOMAIN_data.h5'.
          (e.g. ats_vis_surface_data.h5).
        mesh_filename : str, optional
          Filename for the h5 mesh file.  Default is 'ats_vis_DOMAIN_mesh.h5'.

        Returns
        -------
        self : VisFile object
        """
        self.directory = directory
        self.domain = domain

        self.filename = filename
        if self.filename is None:
            if self.domain is None:
                self.filename = 'ats_vis_data.h5'
            else:
                self.filename = 'ats_vis_{}_data.h5'.format(self.domain)

        self.mesh_filename = mesh_filename
        if self.mesh_filename is None:
            if self.domain is None:
                self.mesh_filename = 'ats_vis_mesh.h5'
            else:
                self.mesh_filename = 'ats_vis_{}_mesh.h5'.format(self.domain)

        if time_unit == 'yr':
            time_factor = 1.0
        elif time_unit == 'noleap':
            time_factor = 365.25 / 365
        elif time_unit == 'd':
            time_factor = 365.25
        elif time_unit == 'hr':
            time_factor = 365.25 * 24
        elif time_unit == 's':
            time_factor = 365.25 * 24 * 3600
        else:
            raise ValueError("Invalid time unit '{}': must be one of 'yr', 'noleap', 'd', 'hr', or 's'".format(time_unit))
        self.time_factor = time_factor
        self.time_unit = time_unit

        self.fname = os.path.join(self.directory, self.filename)
        if not os.path.isfile(self.fname):
            raise RuntimeError("Cannot load ATS XDMF h5 file at: {}".format(self.fname))
        self.d = h5py.File(self.fname,'r')
        self.loadTimes()
        self.map = None
        
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.d.close()

    def loadTimes(self):
        """(Re-)loads the list of cycles and times."""
        a_field = next(iter(self.d.keys()))
        self.cycles = list(sorted(self.d[a_field].keys(), key=int))
        self.times = np.array([self.d[a_field][cycle].attrs['Time'] for cycle in self.cycles]) * self.time_factor

    def filterIndices(self, indices):
        """Filter based on the index into the current set of cycles.

        Note that filters are applied sequentially, but can be undone by
        calling load_times().

        Parameters
        ----------
        indices : one of:
          * int : limits to the ith cycle.
          * list(int) : a list of specific indices
          * slice object : slice the cycle list
        """
        assert(len(self.cycles) == len(self.times))

        if type(indices) is int:
            assert(indices < len(self.cycles))
            self.cycles = [self.cycles[indices],]
            self.times = np.array([self.times[indices],])
        elif type(indices) is slice:
            self.cycles = self.cycles[indices]
            self.times = self.times[indices]
        else:
            inds = list(indices)
            assert(max(inds) < len(self.cycles))
            self.cycles = [self.cycles[i] for i in inds]
            self.times = np.array([self.times[i] for i in inds])
        
    def filterCycles(self, cycles):
        """Filter the vis file based on cycles.

        Note that filters are applied sequentially, but can be undone by
        calling loadTimes().

        Parameters
        ----------
        cycles :
          One of:
          * int : limits to one specific cycle
          * list(int) : a list of specific cycles
        """
        raise RuntimeError
        if type(cycles) is int or type(cycles) is str:
            cycles = [cycles,]
        cycles = [str(c) for c in cycles]

        # note this would be faster with np.isin, but we care about order and
        # repetition of cycles here.
        inds = [np.argwhere(self.cycles == c)[0] for c in cycles]
        self.filterIndices(inds)

    def filterTimes(self, times, eps=1.0):
        """Filter the vis file based on cycles.

        Note that filters are applied sequentially, but can be undone by
        calling load_times().

        Parameters
        ----------
        cycles :
          One of:
          * int : limits to one specific cycle, or the last cycle if -1.
          * list(int) : a list of specific cycles
          * slice object : slice the cycle list
        times : optional
          One of:
          * float : a specific time (within eps in seconds).
          * list(float) : a list of specific times (within eps in seconds).
        eps : float
          Tolerance for defining times, in seconds.  Default is 1.
        """
        if type(times) is float:
            times = [times,]

        # note this would be faster with np.isin, but we care about order and
        # repetition of cycles here.
        inds = [np.argwhere(np.isclose(self.times, t, eps))[0] for t in times]
        self.filterIndices(inds)

    def variable(self, vname):
        """Forms a variable name.

        Parameters
        ----------
        vname : str
          Base variable name, e.g. 'pressure'

        Returns
        -------
        variable_name : str
          Variable name mangled like it is used in Amanzi/ATS.  Something like
          'DOMAIN-vname.cell.0'
        """
        if self.domain and '-' not in vname:
            vname = self.domain + '-' + vname
        if '.' not in vname:
            vname = vname + '.cell.0'
        return vname

    def _get(self, vname, cycle):
        """Private get: assumes vname is fully resolved, and does not deal with maps."""
        return self.d[vname][cycle][:,0]
    
    def get(self, vname, cycle):
        """Access a data member.

        Parameters
        ----------
        vname : str
          Base variable name, e.g. 'pressure'
        cycle : int
          Cycle to access.
        
        Returns
        -------
        value : np.array
          Array of values.

        """
        val = self._get(self.variable(vname), cycle)
        if self.map is None:
            return val
        else:
            return reorder(val, self.map)
        return 

    def getArray(self, vname):
        """Access an array of all cycle values.

        Parameters
        ----------
        vname : str
          Base variable name, e.g. 'pressure'

        Returns
        -------
        value : np.ndarray
          Array of values of shape (n_cycles, n_elems)
        """
        vname = self.variable(vname)
        val = np.array([self._get(vname, k) for k in self.cycles])
        if self.map is None:
            return val
        else:
            return reorder(val, self.map)
            
    
    def loadMesh(self, cycle=None, order=None, shape=None, columnar=False, round=5):
        """Load and reorder centroids and volumes of mesh.

        Parameters
        ----------
        cycle : int, optional
          If the mesh deforms, the centroids may change.  If not provided, gives
          the first cycle's value.
        order : list(str), optional
          See arguments to structuredOrdering().  If provided, this
          reorders the data, and all future get() and getArray() calls will
          return data in this order.
        round : int
          Decimal places to round centroids to.  Supports sorting.
        """

        if cycle is None:
            cycle = self.cycles[0]

        centroids = meshElemCentroids(self.directory, self.mesh_filename, cycle, round)
        if order is None and shape is None and not columnar:
            self.map = None
            self.centroids = centroids
       
        else:
            self.centroids, self.map = structuredOrdering(centroids, order, shape, columnar)

        self.volume = self.get('cell_volume', cycle)

    def loadMeshPolygons(self, cycle=None):
        """Load a mesh into 2D polygons."""
        if cycle is None:
            cycle = self.cycles[0]

        self.loadMesh()
        mesh_elems = meshXYZ(self.directory, self.mesh_filename, cycle)
        self.mesh_elem_info = mesh_elems
        self.polygon_coordinates = meshElemPolygons(*mesh_elems)

    def getMeshPolygons(self, edgecolor='k', cmap='jet', linewidth=1):
        polygons = matplotlib.collections.PolyCollection(self.polygon_coordinates, edgecolor=edgecolor, cmap=cmap, linewidths=linewidth)
        return polygons
    

elem_type = {5:'QUAD',
             8:'PRISM',
             9:'HEX',
             4:'TRIANGLE'
             }


def meshXYZ(directory=".", filename="ats_vis_mesh.h5", key=None):
    """Reads a mesh nodal coordinates and connectivity.

    Note this only currently works for fixed structure meshes, i.e. not
    arbitrary polyhedra.

    Parameters
    ----------
    directory : str, optional
      Directory to read mesh files from.  Default is '.'
    filename : str, optional
      Mesh filename. Default is the Amanzi/ATS default name, 'ats_vis_mesh.h5'
    key : str, optional
      Key of mesh within the file.  This is the cycle number, defaults to the
      first mesh found in the file.

    Returns
    -------
    elemtype : str
      One of 'QUAD', 'PRISM', 'HEX', or 'TRIANGLE' if typed mesh, or 'MIXED' if mesh has more than one types including 'NSIDED' and 'NFACED' 
    coords : np.ndarray
      2D nodal coordinate array.  Shape is (n_nodes, dimension).
    conn : np.ndarray
      2D connection array.  Shape is (n_elem, n_nodes_per_elem + 1), where the
      0th entry in each row is the element type enum, and the remainder of the
      entries are the indices into the nodal array.

    """
    with h5py.File(os.path.join(directory, filename), 'r') as dat:
        if key is None:
            key = next(iter(dat.keys()))

        mesh = dat[key]['Mesh']
        elem_conn = mesh['MixedElements'][:,0]
        coords = coords = mesh['Nodes'][:]
        #coords = dict(zip(mesh['NodeMap'][:,0], mesh['Nodes'][:]))
        elem_type, conns = read_conn(elem_conn)

    return elem_type, coords, conns


def meshElemPolygons(etype, coords, conn):
    """Given mesh info that is a bunch of HEXes, make polygons for 2D plotting."""
    if etype != 'HEX':
        raise RuntimeError("Only works for Hexs")

    y_mean = np.array([c[1] for c in coords.values()]).mean()
    
    coords2 = np.array([[coords[i][0::2] for i in c[1:] if coords[i][1] > y_mean] for c in conn])
    try:
        assert coords2.shape[2] == 2
        assert coords2.shape[1] == 4
    except AssertionError:
        print(coords2.shape)
        for c in conn:
            if len(c) != 9:
                print(c)
                raise RuntimeError("what is a conn?")
            coords3 = np.array([coords[i][:] for i in c[1:] if coords[i][1] > y_mean])
            if coords3.shape[0] != 4:
                print(coords)
                raise RuntimeError("Unable to squash to 2D")

    # reorder anti-clockwise
    for i,c in enumerate(coords2):
        centroid = c.mean(axis=0)
        def angle(p1):
            a1 = np.arctan2((p1[1]-centroid[1]),(p1[0]-centroid[0]))
            return a1

        c2 = np.array(sorted(c,key=angle))
        coords2[i] = c2

    return coords2

def meshElemCentroids(directory=".", filename="ats_vis_mesh.h5", key=None, round=5):
    """Reads and calculates mesh element centroids.

    Note this only currently works for fixed structure meshes, i.e. not
    arbitrary polyhedra.

    Parameters
    ----------
    directory : str, optional
      Directory to read mesh files from.  Default is '.'
    filename : str, optional
      Mesh filename. Default is the Amanzi/ATS default name, 'ats_vis_mesh.h5'
    key : str, optional
      Key of mesh within the file.  This is the cycle number, defaults to the
      first mesh found in the file.
    round : int, optional
      Number of decimals to round to -- this avoids roundoff issues in sorting.
      Default is 5.

    Returns
    -------
    centroids : np.ndarray
      2D nodal coordinate array.  Shape is (n_elems, dimension).

    """
    elem_type, coords, conn = meshXYZ(directory, filename, key)
    centroids = np.zeros((len(conn),3),'d')
    for i,elem in enumerate(conn):
        elem_coords = np.array([coords[gid] for gid in elem[1:]])
        elem_z = np.mean(elem_coords, axis=0)
        centroids[i,:] = elem_z
    return np.round(centroids, round)
    

def structuredOrdering(coordinates, order=None, shape=None, columnar=False):
    """Reorders coordinates in a natural ordering for structured meshes.

    This generates a mapping from an unsorted, unraveled list of
    elements to either a (if shape is None) sorted, unraveled list of
    elements or a (if shape is provided) len(order)+1 dimensional
    array of elements that can be used for structured plotting,
    sorting output, post-processing, and more.

    Parameters
    ----------
    coordinates : np.ndarray
      The 2D array of coordinates, shape (n_coordinates, dimension).
    order : list
      An ordering given to sort(), where headings are 'x', 'y', and potentially
      'z' for 3D meshes.  Note omitted headings always go first.  See below for
      common examples.  May be omitted if columnar == True
    shape : list, optional
      If provided, a list of the same length as order, providing the
      number of elements in that dimension.
    columnar : bool, optional
      If True, this sets order = ['x', 'y', 'z'] and shape is guessed by
      assuming (x,y) coordinates are constant in map-view and cells
      are vertical columns, a common Amanzi/ATS mesh layout.

    Returns
    -------
    ordered_coordinates : np.ndarray
      The re-ordered coordinates, shape (n_coordinates, dimension).
    map : np.ndarray(int)
      Indices of the new coordinates in the old array.  If shape is
      not provided, this is a 1D array and 
        ordered_coordinates[i] == coordinates[map[i]]
      If shape is provided or guess_shape is True, this is a
      len(shape)+1-D array and: 
        ordered_coordinates[i,...,k] == coordinates[map[i,...,k]]

    Examples
    --------

    Sort a column of 100 unordered cells into a 1D sorted array.  The
    input and output are both of shape (100,3).

      > ordered_centroids, map = structuredOrdering(centroids, list())

    Sort a logically structured transect of size NX x NY x NZ =
    (100,1,20), where x is structured and z may vary as a function of
    x.  Both input and output are of shape (2000, 3), but the output
    is sorted with each column appearing sequentially and the
    z-dimension fastest-varying.  map is of shape (2000,).
    
      > ordered_centroids, map = structuredOrdering(centroids, ['z',])

    Do the same, but this time reshape into a 2D array.  Now the
    ordered_centroids are of shape (100, 20, 3), and the map is of
    shape (100, 20).
    
      > ordered_centroids, map = structuredOrdering(centroids, ['z',], [20,])

    Do the same as above, but detect the shape.  This works only
    because the mesh is columnar.

      > ordered_centroids, map = structuredOrdering(centroids, columnar=True)

    Sort a 3D map-view "structured-in-z" mesh into arbitrarily-ordered
    x and y columns.  Assume there are 1000 map-view triangles, each
    extruded 20 cells deep.  The input is is of shape (20000, 3) and
    the output is of shape (1000, 20, 3).

      > ordered_centroids, map = structuredOrdering(centroids, columnar=True)

    Note that map can be used with the reorder() function to place
    data in this ordering.

    """
    if columnar:
        order = ['x', 'y', 'z',]
        
    
    # Surely there is a cleaner way to do this in numpy?
    # The current approach packs, sorts, and unpacks.
    if (coordinates.shape[1] == 3):
        coords_a = np.array([(i,coordinates[i,0],coordinates[i,1],coordinates[i,2])
                             for i in range(coordinates.shape[0])],
                            dtype=[('id',int),('x',float),('y',float),('z',float)])
    elif (coordinates.shape[1] == 2):
        coords_a = np.array([(i,coordinates[i,0],coordinates[i,1])
                             for i in range(coordinates.shape[0])],
                            dtype=[('id',int),('x',float),('y',float)])
    coords_a.sort(order=order)
    map = coords_a['id']
    if (coordinates.shape[1] == 3):
        ordered_coordinates = np.array([coords_a['x'], coords_a['y'], coords_a['z']]).transpose()
    else:
        ordered_coordinates = np.array([coords_a['x'], coords_a['y']]).transpose()
        
    if columnar:
        # try to guess the shape based on new-found contiguity
        n_cells_in_column = 0
        xy = ordered_coordinates[0,0:2]
        while n_cells_in_column < ordered_coordinates.shape[0] and \
              np.allclose(xy, ordered_coordinates[n_cells_in_column,0:2], 0., 1.e-5):
            n_cells_in_column += 1
        shape = [n_cells_in_column,]
     

    if shape is not None:
        new_shape = (-1,) + tuple(shape)
        coord_shape = new_shape+(3,)
        ordered_coordinates = np.reshape(ordered_coordinates, coord_shape)
        map = np.reshape(map, new_shape)

        if map.shape[0] == 1:
            map = map[0]
            ordered_coordinates = ordered_coordinates[0]

    return ordered_coordinates, map


def reorder(data, map):
    """Re-orders values and arrays according to a mesh reordering.

    Parameters
    ----------
    data : np.ndarray
      The data, i.e. provided by VisFile.get() or VisFile.getArray() 
    map : np.ndarray
      A reordering of indices to remap the data based on a map
      returned by structuredOrdering()

    Returns
    -------
    data : np.ndarray
      The re-ordered data.

    """
    flatten = (len(data.shape) == 1)
    if flatten:
        data = np.expand_dims(data, 0)

    # unravel the map, reorder, then reshape back into map's shape
    map_shape = map.shape
    data = data[:, map.ravel()].reshape((-1,)+map_shape)

    # if one cycle, get one data
    if flatten:
        data = data[0]

    return data


elem_type = {3:'POLYGON',
             5:'QUAD',
             8:'PRISM',
             9:'HEX',
             4:'TRIANGLE',
             16:'POLYHEDRON'
             }

elem_typed_node_counts = { 'QUAD' : 4,
                     'PRISM' : 6,
                     'HEX' : 8,
                     'TRIANGLE' : 3
                     }


def read_conn(elem_conn):
    """Reads an array, called MixedElements in the HDF5 file, to get conn"""
    i = 0
    etypes = []
    conns = []
    while i < len(elem_conn):
        etype, conn, i = read_element_dirty(elem_conn,i)
        etypes.append(etype)
        conns.append(conn)
    if len(set(etypes))==1:
        elem_type = set(etypes)[0]
    else: 
        elem_type = 'MIXED'
    return elem_type, conns


def read_element_dirty(elem_conn, i):
    """Reads the element at location i,
    
    returns etype, nodeids, new_i

    Note this is called dirty because it does not _properly_ deal with
    NFACED objects, but instead just returns a set of unique nodes
    that are in the element (i.e. it has no concept of faces).

    """
    try:
        etype = elem_type[elem_conn[i]]
    except KeyError:
        raise RuntimeError(f'This reader is not implemented for elements of type {elem_conn[i]} -- what type is this?')
    if etype == 'POLYGON':
        return 'POLYGON', *read_polygon_element(elem_conn, i+1)
    elif etype == 'POLYHEDRON':
        return 'POLYHEDRON', *read_polyhedron_element_dirty(elem_conn, i+1)
    else:
        return etype, *read_typed_element(elem_conn, etype, i+1)


def read_polygon_element(elem_conn, i):
    n_nodes = elem_conn[i]
    nodes = elem_conn[i+1:i+1+n_nodes]
    return nodes, i+1+n_nodes


def read_typed_element(elem_conn, etype, i):
    n_nodes = elem_typed_node_counts[etype]
    nodes = elem_conn[i:i+n_nodes]
    return nodes, i+n_nodes


def read_polyhedron_element_dirty(elem_conn, i):
    n_faces = elem_conn[i]
    i = i + 1
    elem_nodes = set()
    for j in range(n_faces):
        (fnodes, i) = read_polygon_element(elem_conn, i)
        elem_nodes = elem_nodes.union(fnodes)
    return list(elem_nodes), i
    
        
    

