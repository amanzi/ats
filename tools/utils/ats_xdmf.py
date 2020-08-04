"""Functions for parsing Amanzi/ATS XDMF visualization files."""
import sys,os
import numpy as np
import h5py

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
          Filename of h5 vis file.  Default is 'visdump_DOMAIN_data.h5'.
          (e.g. visdump_surface_data.h5).
        mesh_filename : str, optional
          Filename for the h5 mesh file.  Default is 'visdump_DOMAIN_mesh.h5'.

        Returns
        -------
        self : VisFile object
        """
        self.directory = directory
        self.domain = domain

        self.filename = filename
        if self.filename is None:
            if self.domain is None:
                self.filename = 'visdump_data.h5'
            else:
                self.filename = 'visdump_{}_data.h5'.format(self.domain)

        self.mesh_filename = mesh_filename
        if self.mesh_filename is None:
            if self.domain is None:
                self.mesh_filename = 'visdump_mesh.h5'
            else:
                self.mesh_filename = 'visdump_{}_mesh.h5'.format(self.domain)

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
        
        self.d = h5py.File(os.path.join(self.directory, self.filename))
        self.loadTimes()
        self.map = None
        
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.d.close()

    def loadTimes(self):
        """(Re-)loads the list of cycles and times."""
        a_field = next(iter(self.d.keys()))
        self.cycles = np.array(sorted(self.d[a_field].keys(), key=int))
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
        self.cycles = self.cycles[indices]
        self.times = self.times[indices]

        
    def filterCycles(self, cycles):
        """Filter the vis file based on cycles.

        Note that filters are applied sequentially, but can be undone by
        calling loadTimes().

        Parameters
        ----------
        cycles :
          One of:
          * int : limits to one specific cycle, or the last cycle if -1.
          * list(int) : a list of specific cycles
        """
        if type(cycles) is int:
            cycles = [cycles,]

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
          See arguments to mesh.structuredOrdering().  If provided, this
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

            
    
elem_type = {5:'QUAD',
             8:'PRISM',
             9:'HEX',
             4:'TRIANGLE'
             }

def meshXYZ(directory=".", filename="visdump_mesh.h5", key=None):
    """Reads a mesh nodal coordinates and connectivity.

    Note this only currently works for fixed structure meshes, i.e. not
    arbitrary polyhedra.

    Parameters
    ----------
    directory : str, optional
      Directory to read mesh files from.  Default is '.'
    filename : str, optional
      Mesh filename. Default is the Amanzi/ATS default name, 'visdump_mesh.h5'
    key : str, optional
      Key of mesh within the file.  This is the cycle number, defaults to the
      first mesh found in the file.

    Returns
    -------
    etype : str
      One of 'QUAD', 'PRISM', 'HEX', or 'TRIANGLE'.  Note 'NSIDED' and 'NFACED' 
      are not yet supported.
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

        etype = elem_type[elem_conn[0]]
        if (etype == 'PRISM'):
            nnodes_per_elem = 6
        elif (etype == 'HEX'):
            nnodes_per_elem = 8
        elif (etype == 'QUAD'):
            nnodes_per_elem = 4
        elif (etype == 'TRIANGLE'):
            nnodes_per_elem = 3

        if len(elem_conn) % (nnodes_per_elem + 1) != 0:
            raise ValueError('This reader only processes single-element-type meshes.')
        n_elems = int(len(elem_conn) / (nnodes_per_elem+1))
        coords = dict(zip(mesh['NodeMap'][:,0], mesh['Nodes'][:]))

    conn = elem_conn.reshape((n_elems, nnodes_per_elem+1))
    if (np.any(conn[:,0] != elem_conn[0])):
        raise ValueError('This reader only processes single-element-type meshes.')
    return etype, coords, conn


def meshElemCentroids(directory=".", filename="visdump_mesh.h5", key=None, round=5):
    """Reads and calculates mesh element centroids.

    Note this only currently works for fixed structure meshes, i.e. not
    arbitrary polyhedra.

    Parameters
    ----------
    directory : str, optional
      Directory to read mesh files from.  Default is '.'
    filename : str, optional
      Mesh filename. Default is the Amanzi/ATS default name, 'visdump_mesh.h5'
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
    etype, coords, conn = meshXYZ(directory, filename, key)

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
    # one cycle or many?
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
