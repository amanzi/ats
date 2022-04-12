ATS Native XML Input Specification V-dev
****************************************


.. contents:: **Table of Contents**
   :depth: 2

  
Syntax of the Specification
#######################################

* Input specification for each ParameterList entry consists of two parts.  
  First, a bulleted list defines the usage syntax and available options.  
  This is followed by example snipets of XML code to demonstrate usage.

* In many cases, the input specifies data for a particular parameterized model, and ATS 
  supports a number of parameterizations.  
  For example, initial data might be uniform (the value is required), or linear in y (the value 
  and its gradient are required).  
  Where ATS supports a number of parameterized models for quantity Z, the available 
  models will be listed by name, and then will be described in the subsequent section.  
  For example, the specification for an `"X`" list might begin with the following:

  * `"Y`" ``[string]`` **"default_value"**, `"other`", `"valid`", `"options`"

  * Z ``[Z-spec]`` Model for Z, choose exactly one of the following: (1) `"z1`", or (2) `"z2`" (see below) 

Here, an `"X`" is defined by a `"Y`" and a `"Z`".  
The `"Y`" is a string parameter but the `"Z`" is given by a model (which will require its own set of parameters).
The options for `"Z`" will then be described as a spec:

 * `"z1`" applies model z1.  Requires `"z1a`" ``[string]``

 * `"z2`" applies model z2.  Requires `"z2a`" ``[double]`` and `"z2b`" ``[int]``

An example of using such a specification:

.. code-block:: xml

    <ParameterList name="X">
      <Parameter name="Y" type="string" value="hello"/>
      <ParameterList name="z2">
        <Parameter name="z2a" type="double" value="0.7"/>
        <Parameter name="z2b" type="int" value="3"/>
      </ParameterList>   
    </ParameterList>   
 
Here, the user is defining X with Y="hello", and Z will be a z2 constructed with z2a=0.7 and z2b=3.

Conventions:

* Reserved keywords and labels are `"quoted and italicized`" -- these
  labels or values of parameters in user-generated input files must
  match (using XML matching rules) the specified or allowable values.

* User-defined labels are indicated with ALL-CAPS, and are meant to
  represent a typical name given by a user - these can be names or
  numbers or whatever serves best the organization of the user input
  data.

* Bold values are default values, and are used if the Parameter
  is not provided.


Symbol Index
############

.. include:: symbol_table.rst
  
Main
####


ATS's top-level main accepts an XML list including a few required elements.

.. _main-spec:
.. admonition:: main-spec

    * `"mesh`" ``[mesh-typed-spec-list]`` A list of Mesh_ objects.
    * `"regions`" ``[region-typedsublist-spec-list]`` A list of Region_ objects.
    * `"cycle driver`" ``[coordinator-spec]``  See Coordinator_.
    * `"visualization`" ``[visualization-spec-list]`` A list of Visualization_ objects.
    * `"observations`" ``[observation-spec-list]`` An list of Observation_ objects.
    * `"checkpoint`" ``[checkpoint-spec]`` See Checkpoint_.
    * `"PKs`" ``[pk-typedinline-spec-list]`` A list of PK_ objects.
    * `"state`" ``[state-spec]`` See State_.

 

  

Mesh
####
 A list of mesh objects and their domain names.

All processes are simulated on a domain, which is discretized through a mesh.

Multiple domains and therefore meshes can be used in a single simulation, and
multiple meshes can be constructed on the fly.  The top level `"mesh`" is a
list of ``[mesh-typed-spec]`` sublists whose name indicate the mesh or domain
name.

Included in that list is at least one mesh: the `"domain`" mesh.  The
`"domain`" mesh represents the primary domain of simulation -- usually the
subsurface.  Simple, structured meshes may be generated on the fly, or complex
unstructured meshes are provided as Exodus II files.  The `"domain`" mesh list
includes either a `Generated Mesh`_, `Read Mesh File`_, or `Logical Mesh`_ spec, as
described below.

Additionally, a `Surface Mesh`_ may be formed by lifting the surface of a
provided mesh and then flattening that mesh to a 2D surface.  `Column Meshes`_
which split a base mesh into vertical columns of cells for use in 1D models
may also be generated automatically.

Finally, mesh generation is hard and error-prone.  A mesh audit is provided,
which checks for many common geometric and topologic errors in mesh
generation.  This is reasonably fast, even for big meshes, and can be done
through providing a "verify mesh" option.

.. _mesh-typed-spec:
.. admonition:: mesh-typed-spec

    * `"mesh type`" ``[string]`` One of:

      - `"generate mesh`" See `Generated Mesh`_.
      - `"read mesh file`" See `Read Mesh File`_.
      - `"logical`" See `Logical Mesh`_.
      - `"surface`" See `Surface Mesh`_.
      - `"subgrid`" See `Subgrid Meshes`_.
      - `"column`" See `Column Meshes`_.

    * `"_mesh_type_ parameters`" ``[_mesh_type_-spec]`` List of parameters
      associated with the type.
    * `"verify mesh`" ``[bool]`` **false** Perform a mesh audit.
    * `"deformable mesh`" ``[bool]`` **false** Will this mesh be deformed?

    * `"build columns from set`" ``[string]`` **optional** If provided, build
       columnar structures from the provided set.

    * `"partitioner`" ``[string]`` **zoltan_rcb** Method to partition the
      mesh.  Note this only makes sense on the domain mesh.  One of:

      - `"zoltan_rcb`" a "map view" partitioning that keeps columns of cells together
      - `"metis`" uses the METIS graph partitioner
      - `"zoltan`" uses the default Zoltan graph-based partitioner.


Generated Mesh
==============

Generated mesh are by definition structured, with uniform dx, dy, and dz.
Such a mesh is specified by a bounding box high and low coordinate, and a list
of number of cells in each direction.

Specified by `"mesh type`" of `"generate mesh`".

.. _mesh-generate-mesh-spec:
.. admonition:: mesh-generate-mesh-spec

    * `"domain low coordinate`" ``[Array(double)]`` Location of low corner of domain
    * `"domain high coordinate`" ``[Array(double)]`` Location of high corner of domain
    * `"number of cells`" ``[Array(int)]`` the number of uniform cells in each coordinate direction

Example:

.. code-block:: xml

   <ParameterList name="mesh">
     <ParameterList name="domain">
       <Parameter name="mesh type" type="string" value="generate mesh"/>
       <ParameterList name="generate mesh parameters"/>
         <Parameter name="number of cells" type="Array(int)" value="{{100, 1, 100}}"/>
         <Parameter name="domain low coordinate" type="Array(double)" value="{{0.0, 0.0, 0.0}}" />
         <Parameter name="domain high coordinate" type="Array(double)" value="{{100.0, 1.0, 10.0}}" />
       </ParameterList>
     </ParameterList>   
   </ParameterList>   


Read Mesh File
==============

Meshes can be pre-generated in a multitude of ways, then written to file, and
loaded in ATS. Note that in the case of an Exodus II mesh file, the suffix of
the serial mesh file must be .exo and the suffix of the parallel mesh file must
be .par.  When running in serial the code will read this the indicated file
directly.  When running in parallel with a prepartitioned mesh, the suffix is
.par and the code will instead read the partitioned files that have been
generated with a Nemesis tool and named as filename.par.N.r where N is the
number of processors and r is the rank.  When running in parallel and the
suffix is .exo, the code will partition automatically the serial file.

Specified by `"mesh type`" of `"read mesh file`".

.. _mesh-read-mesh-file-spec:
.. admonition:: mesh-read-mesh-file-spec

    * `"file`" ``[string]`` filename of a pre-generated mesh file

Example:

.. code-block:: xml

   <ParameterList name="mesh">
     <ParameterList name="domain">
       <Parameter name="mesh type" type="string" value="read mesh file"/>
       <ParameterList name="read mesh file parameters">
         <Parameter name="file" type="string" value="mesh_filename.exo"/>
       </ParameterList>
       <Parameter name="verify mesh" type="bool" value="true" />
     </ParameterList>
   </ParameterList>


Logical Mesh
============

Logical meshes are meshes for whom nodal coordinates may not be
specified, but sufficient information about the geometry of the
conceptual domain can be specified to allow solving problems.  This
allows for the conceptual generation of domains that "act" like a mesh
and can be used like a mesh, but don't fit MSTK's view of an
unstructured mesh.

This is an active research and development area, and is used most
frequently for river networks, root networks, and crack networks.

Specified by `"mesh type`" of `"logical`".

.. todo::
   WIP: add spec!

.. _mesh-logical-spec:
.. admonition:: mesh-logical-spec

    Not yet completed...
   
Surface Mesh
============

To lift a surface off of the mesh, a side-set specifying all surface faces
must be given.  These faces are lifted locally, so the partitioning of the
surface cells will be identical to the partitioning of the subsurface faces
that correspond to these cells.  All communication and ghost cells are set up.
The mesh is flattened, so all surface faces must have non-zero area when
projected in the z-direction.  No checks for holes are performed.  Surface
meshes may similarly be audited to make sure they are reasonable for
computation.

Specified by `"mesh type`" of `"surface`".

.. _mesh-surface-spec:
.. admonition:: mesh-surface-spec

    ONE OF

    * `"surface sideset name`" ``[string]`` The Region_ name containing all surface faces.

    OR

    * `"surface sideset names`" ``[Array(string)]`` A list of Region_ names containing the surface faces.

    END

    * `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
    * `"export mesh to file`" ``[string]`` **optional** Export the lifted
      surface mesh to this filename.

Example:

.. code-block:: xml

    <ParameterList name="mesh" type="ParameterList">
      <ParameterList name="surface" type="ParameterList">
        <Parameter name="mesh type" type="string" value="surface" />
        <ParameterList name="surface parameters" type="ParameterList">
          <Parameter name="surface sideset name" type="string" value="{surface_region}" />
          <Parameter name="verify mesh" type="bool" value="true" />
          <Parameter name="export mesh to file" type="string" value="surface_mesh.exo" />
        </ParameterList>
      </ParameterList>
      <ParameterList name="domain" type="ParameterList">
        <Parameter name="mesh type" type="string" value="read mesh file" />
        <ParameterList name="read mesh file parameters" type="ParameterList">
          <Parameter name="file" type="string" value="../data/open-book-2D.exo" />
          <Parameter name="format" type="string" value="Exodus II" />
        </ParameterList>
      </ParameterList>
    </ParameterList>


Aliased Mesh
============

Aliased domains are simply domains that share a mesh with another domain.  For
instance, one might find it useful to define both a "surface water" and a
"snow" domain that share a common "surface" mesh.  In that case, typically the
"surface" domain would point to the "surface" mesh, and the "snow" domain would
be an "aliased" domain whose target is the "surface" mesh.

Specified by `"mesh type`" of `"aliased`".

.. _mesh-aliased-spec:
.. admonition:: mesh-aliased-spec

    * `"target`" ``[string]`` Mesh that this alias points to.


Subgrid Meshes
==============

A collection of meshes formed by associating a new mesh with each entity of a
region.  Used for a few cases, including generating a 1D column for each
surface face of a semi-structured subsurface mesh, or for hanging logical
meshes off of each surface cell as a subgrid model, etc.

The subgrid meshes are then named `"MESH_NAME_X"` for each X, which is an
entity local ID, in a provided region of the provided entity type.

Specified by `"mesh type`" of `"subgrid`".

.. _mesh-subgrid-spec:
.. admonition:: mesh-subgrid-spec

    * `"subgrid region name`" ``[string]`` Region on which each subgrid mesh will be associated.
    * `"entity kind`" ``[string]`` One of `"cell`", `"face`", etc.  Entity of the
      region (usually `"cell`") on which each subgrid mesh will be associated.
    * `"parent domain`" ``[string]`` **domain** Mesh which includes the above region.
    * `"flyweight mesh`" ``[bool]`` **False** NOT YET SUPPORTED.  Allows a single
      mesh instead of one per entity.

.. todo::
   WIP: Add examples (intermediate scale model, transport subgrid model)


Column Meshes
=============

.. warning::
   Note these are rarely if ever created manually by a user.  Instead use
   `Subgrid Meshes`_, which generate a column mesh spec for every face
   of a set.

Specified by `"mesh type`" of `"column`".

.. _mesh-column-spec:
.. admonition:: mesh-column-spec

    * `"parent domain`" ``[string]`` The name of the 3D mesh from which columns are generated.
      Note that the `"build columns from set`" parameter must be set in that mesh.
    * `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
    * `"deformable mesh`" ``[bool]`` **false**  Used for deformation PKs to allow non-const access.
    * `"entity LID`" ``[int]`` Local ID of the surface cell that is the top of the column.

Example:

.. code-block:: xml

    <ParameterList name="mesh" type="ParameterList">
      <ParameterList name="column" type="ParameterList">
        <ParameterList name="column parameters" type="ParameterList">
          <Parameter name="parent domain" type="string" value="domain" />
          <Parameter name="entity LID" type="int" value="0" />
        </ParameterList>
      </ParameterList>
      <ParameterList name="domain" type="ParameterList">
        <Parameter name="mesh type" type="string" value="read mesh file" />
        <ParameterList name="read mesh file parameters" type="ParameterList">
          <Parameter name="file" type="string" value="../data/open-book-2D.exo" />
          <Parameter name="format" type="string" value="Exodus II" />
        </ParameterList>
      </ParameterList>
    </ParameterList>





Region
######
  A geometric or discrete subdomain of the full domain.

Regions are geometrical constructs used to define subsets of
the computational domain in order to specify the problem to be solved, and the
output desired. Regions may represents zero-, one-, two- or three-dimensional
subsets of physical space.  For a three-dimensional problem, the simulation
domain will be a three-dimensional region bounded by a set of two-dimensional
regions.  If the simulation domain is N-dimensional, the boundary conditions
must be specified over a set of regions are (N-1)-dimensional.

Region specs are **not** denoted by a "type" parameter for legacy reasons.
Instead, they take a single sublist whose name defines the type.

``[region-typedsublist-spec]``


.. warning:: Surface files contain labeled triangulated face sets.  The user is
    responsible for ensuring that the intersections with other surfaces in the
    problem, including the boundaries, are *exact* (*i.e.* that surface
    intersections are *watertight* where applicable), and that the surfaces are
    contained within the computational domain.  If nodes in the surface fall
    outside the domain, the elements they define are ignored.

    Examples of surface files are given in the *Exodus II* file format here.

.. warning:: Region names must NOT be repeated.

Example:

.. code-block:: xml

   <ParameterList>  <!-- parent list -->
     <ParameterList name="regions">
       <ParameterList name="TOP SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 5}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 8}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="MIDDLE SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 3}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 5}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BOTTOM SECTION">
         <ParameterList name="region: box">
           <Parameter name="low coordinate" type="Array(double)" value="{2, 3, 0}"/>
           <Parameter name="high coordinate" type="Array(double)" value="{4, 5, 3}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="INFLOW SURFACE">
         <ParameterList name="region: labeled set">
           <Parameter name="label"  type="string" value="sideset_2"/>
           <Parameter name="file"   type="string" value="F_area_mesh.exo"/>
           <Parameter name="format" type="string" value="Exodus II"/>
           <Parameter name="entity" type="string" value="face"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="OUTFLOW PLANE">
         <ParameterList name="region: plane">
           <Parameter name="point" type="Array(double)" value="{0.5, 0.5, 0.5}"/>
           <Parameter name="normal" type="Array(double)" value="{0, 0, 1}"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="BLOODY SAND">
         <ParameterList name="region: color function">
           <Parameter name="file" type="string" value="F_area_col.txt"/>
           <Parameter name="value" type="int" value="25"/>
         </ParameterList>
       </ParameterList>
       <ParameterList name="FLUX PLANE">
         <ParameterList name="region: polygon">
           <Parameter name="number of points" type="int" value="5"/>
           <Parameter name="points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                                  0.5, -0.5, -0.5,
                                                                  0.8, 0.0, 0.0,
                                                                  0.5,  0.5, 0.5,
                                                                 -0.5, 0.5, 0.5}"/>
          </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>

In this example, *TOP SECTION*, *MIDDLE SECTION* and *BOTTOM SECTION*
are three box-shaped volumetric regions. *INFLOW SURFACE* is a
surface region defined in an Exodus II-formatted labeled set
file and *OUTFLOW PLANE* is a planar region. *BLOODY SAND* is a volumetric
region defined by the value 25 in color function file.




.. contents:: **Region Types**
   :local:


All
===
  A region consisting of all entities on a mesh.

No parameters required.

``[region-all-spec]``

   * `"empty`" ``[bool]`` **True** This is simply here to avoid issues with
       empty lists.  The better solution is to rewrite the Region spec
       completely to make it consistent with all other typed specs in Amanzi.

Example:

.. code-block:: xml

   <ParameterList name="domain">  <!-- parent list -->
     <ParameterList name="region: all">
     </ParameterList>
   </ParameterList>
 



Box
===
 RegionBox: a rectangular region in space, defined by two corners

List *region: box* defines a region bounded by coordinate-aligned
planes. Boxes are allowed to be of zero thickness in only one
direction in which case they are equivalent to planes.

.. _region-box-spec:
.. admonition:: region-box-spec

    * `"low coordinate`" ``[Array(double)]`` Location of the boundary point with the lowest coordinates.
    * `"high coordinate`" ``[Array(double)]`` Location of the boundary points with the highest coordinates.

Example:

.. code-block:: xml

   <ParameterList name="WELL">  <!-- parent list -->
     <ParameterList name="region: box">
       <Parameter name="low coordinate" type="Array(double)" value="{-5.0,-5.0, -5.0}"/>
       <Parameter name="high coordinate" type="Array(double)" value="{5.0, 5.0,  5.0}"/>
     </ParameterList>
   </ParameterList>
  



Plane
=====
 RegionPlane: A planar (infinite) region in space, defined by a point and a normal.

List *region: plane* defines a plane using a point lying on the plane and normal to the plane.

.. _region-plane-spec:
.. admonition:: region-plane-spec

    * `"normal`" ``[Array(double)]`` Normal to the plane.
    * `"point`" ``[Array(double)]`` Point in space.

Example:

.. code-block:: xml

   <ParameterList name="TOP_SECTION"> <!-- parent list -->
     <ParameterList name="region: plane">
       <Parameter name="point" type="Array(double)" value="{2, 3, 5}"/>
       <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-05"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>




Labeled Set
===========
 RegionLabeledSet: A region defined by a set of mesh entities in a mesh file

The list *region: labeled set* defines a named set of mesh entities
existing in an input mesh file. This is the same file that contains
the computational mesh. The name of the entity set is given
by *label*.  For example, a mesh file in the Exodus II
format can be processed to tag cells, faces and/or nodes with
specific labels, using a variety of external tools. Regions based
on such sets are assigned a user-defined label for Amanzi, which may
or may not correspond to the original label in the exodus file.
Note that the file used to express this labeled set may be in any
Amanzi-supported mesh format (the mesh format is specified in the
parameters for this option).  The *entity* parameter may be
necessary to specify a unique set.  For example, an Exodus file
requires *cell*, *face* or *node* as well as a label (which is
an integer).  The resulting region will have the dimensionality 
associated with the entities in the indicated set. 

.. _region-labeled-set-spec:
.. admonition:: region-labeled-set-spec

    * `"label`" ``[string]`` Set per label defined in the mesh file.
    * `"file`" ``[string]`` File name.
    * `"format`" ``[string]`` Currently, we only support mesh files in the "Exodus II" format.
    * `"entity`" ``[string]`` Type of the mesh object (cell, face, etc).

Example:

.. code-block:: xml

   <ParameterList name="AQUIFER">
     <ParameterList name="region: labeled set">
       <Parameter name="entity" type="string" value="cell"/>
       <Parameter name="file" type="string" value="porflow4_4.exo"/>
       <Parameter name="format" type="string" value="Exodus II"/>
       <Parameter name="label" type="string" value="1"/>
     </ParameterList>
   </ParameterList>




Function Color
==============
 RegionFunctionColor: A region defined by the value of an indicator function in a file.

The list *region: color function* defines a region based a specified integer
color, *value*, in a structured color function file, *file*.  The format of
the color function file is given below in the "Tabulated function file format"
section. As shown in the file, the color values may be specified at the nodes
or cells of the color function grid. A computational cell is assigned the
'color' of the data grid cell containing its cell centroid (cell-based colors)
or the data grid nearest its cell-centroid (node-based colors). Computational
cells sets are then built from all cells with the specified color *Value*.

In order to avoid, gaps and overlaps in specifying materials, it is strongly
recommended that regions be defined using a single color function file.

.. _region-color-function-spec:
.. admonition:: region-color-function-spec

    * `"file`" ``[string]`` File name containing color function.
    * `"value`" ``[int]`` Color that defines the set in the tabulated function file.

Example:

.. code-block:: xml

   <ParameterList name="SOIL_TOP">
     <ParameterList name="region: color function">
       <Parameter name="file" type="string" value="geology_resamp_2D.tf3"/>
       <Parameter name="value" type="int" value="1"/>
     </ParameterList>
   </ParameterList>




Point
=====
 RegionPoint: a point in space.

List *region: point* defines a point in space. 
This region consists of cells containing this point.

.. _region-point-spec:
.. admonition:: region-point-spec

    * `"coordinate`" ``[Array(double)]`` Location of point in space.

Example:

.. code-block:: xml

   <ParameterList name="DOWN_WIND150"> <!-- parent list defining the name -->
     <ParameterList name="region: point">
       <Parameter name="coordinate" type="Array(double)" value="{-150.0, 0.0, 0.0}"/>
     </ParameterList>
   </ParameterList>




Logical
=======
 RegionLogical: A region defined by a logical operation on one or two other regions

The list *region: logical* defines logical operations on regions allow for
more advanced region definitions. At this time the logical region allows for
logical operations on a list of regions.  *union* and *intersection* are
self-evident. In the case of *subtraction*, subtraction is performed from the
first region in the list.  The *complement* is a special case in that it is
the only case that operates on single region, and returns the complement to it
within the domain ENTIRE_DOMAIN.  Currently, multi-region booleans are not
supported in the same expression.

.. _region-logical-spec:
.. admonition:: region-logical-spec

    * `"operation`" ``[string]`` defines operation on the list of regions.
      One of: `"union`", `"intersect`", `"subtract`", `"complement`"
    * `"regions`" ``[Array(string)]`` specifies the list of involved regions.

Example:

.. code-block:: xml

  <ParameterList name="LOWER_LAYERs">
    <ParameterList name="region: logical">
      <Parameter name="operation" type="string" value="union"/>
      <Parameter name="regions" type="Array(string)" value="{Middle1, Middle2, Bottom}"/>
    </ParameterList>
  </ParameterList>




Polygon
=======
 RegionPolygon: A closed polygonal segment of a plane.

The list *region: polygon* defines a polygonal region on which mesh faces and
nodes can be queried. NOTE that one cannot ask for cells in a polygonal surface
region. In 2D, the polygonal region is a line and is specified by 2 points.
In 3D, the polygonal region is specified by an arbitrary number of points.
In both cases the point coordinates are given as a linear array. The polygon
can be non-convex.

This provides a set of faces with a normal for computing flux.

The polygonal surface region can be queried for a normal. In 2D, the normal is
defined as [Vy,-Vx] where [Vx,Vy] is the vector from point 1 to point 2.
In 3D, the normal of the polygon is defined by the order in which points 
are specified.

``[region-polygon-spec]``

* `"number of points`" ``[int]`` Number of polygon points.

* `"points`" ``[Array(double)]`` Point coordinates in a linear array. 

Example:

.. code-block:: xml

   <ParameterList name="XY_PENTAGON">
     <ParameterList name="region: polygon">
       <Parameter name="number of points" type="int" value="5"/>
       <Parameter name="points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                              0.5, -0.5, -0.5,
                                                              0.8, 0.0, 0.0,
                                                              0.5,  0.5, 0.5,
                                                             -0.5, 0.5, 0.5}"/>
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-3"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>




Enumerated
==========
 RegionEnumerated: A region enumerated as a list of IDs.

List *region: enumerated set* defines a set of mesh entities via the list 
of input global ids. Note that global ids are not defined correctly when
parallel mesh is created on a fly.

.. _region-enumerated-spec:
.. admonition:: region-enumerated-spec

    * `"entity`" ``[string]`` Type of the mesh object.  One of: `"cell`", `"face`", `"edge`", `"node`"
    * `"entity gids`" ``[Array(int)]`` List of the global IDs of the entities.
  

Example:

.. code-block:: xml

   <ParameterList name="WELL"> <!-- parent list -->
     <ParameterList name="region: enumerated set">
       <Parameter name="entity" type="string" value="face"/>
       <Parameter name="entity gids" type="Array(int)" value="{1, 12, 23, 34}"/>
     </ParameterList>
   </ParameterList>




Boundary
========
 RegionBoundary:  A region consisting of all entities on the domain boundary

List *region: boundary* defines a set of all boundary faces. 
Using this definition, faces located on the domain boundary are extracted.

.. _region-boundary-spec:
.. admonition:: region-boundary-spec

    * `"entity`" ``[string]`` Type of the mesh object.  Unclear whether this is
      used or can be other things than `"face`"?

Example:

.. code-block:: xml

   <ParameterList name="DOMAIN_BOUNDARY"> <!-- parent list names the region -->
     <ParameterList name="region: boundary">
       <Parameter name="entity" type="string" value="face"/>
     </ParameterList>
   </ParameterList>




Box Volume Fractions
====================
 RegionBoxVolumeFractions: A rectangular region in space, defined by two corner points and normals to sides.

List *region: box volume fraction* defines a region bounded by a box *not* 
aligned with coordinate axes. 
Boxes are allowed to be of zero thickness in only one direction in which case 
they are equivalent to rectangles on a plane or segments on a line.

.. _region-box-volume-fractions-spec:
.. admonition:: region-box-volume-fractions-spec

    * `"corner coordinate`" ``[Array(double)]`` Location of one box corner.
    * `"opposite corner coordinate`" ``[Array(double)]`` Location of the opposite box corner.
    * `"normals`" ``[Array(double)]`` Normals to sides in a linear array. Default is columns of
      the identity matrix. The normals may be scaled arbitrarily but must be orthogonal to
      one another and form the right coordinate frame.

Example:

.. code-block:: xml

   <ParameterList name="BASIN">  <!-- parent list -->
     <ParameterList name="region: box volume fractions">
       <Parameter name="corner coordinate" type="Array(double)" value="{-1.0,-1.0, 1.0}"/>
       <Parameter name="opposite corner coordinate" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
       <Parameter name="normals" type="Array(double)" value="{1.0, 0.0, 0.0
                                                              0.0, 2.0, 0.0,
                                                              0.0, 0.0, 3.0}"/>
     </ParameterList>
   </ParameterList>

This example defines a degenerate box, a square on a surface *z=1*.

 


Line Segment
============
 RegionLineSegment: A line segment, defined by two points in space.

List *region: line segment* desribes a region defined by a line
segment. This region is a set of cells which intersect with a line
segment.  The line segment is allowed to intersect with one or more cells. Zero length
line segments are allowed. The line segment is defined by its ends
points.

.. _region-line-segment-spec:
.. admonition:: region-line-segment-spec

    * `"end coordinate`" ``[Array(double)]`` Location of one end of a line
      segment.
    * `"opposite end coordinate`" ``[Array(double)]`` Location of the opposite
      end of a line segment.

Example:

.. code-block:: xml

   <ParameterList name="WELL"> <!-- parent list -->
      <ParameterList name="region: line segment">
        <Parameter name="end coordinate" type="Array(double)" value="{497542.44, 5393755.77, 0.0}"/>
        <Parameter name="opposite end coordinate" type="Array(double)" value="{497542.44, 5393755.77, 100.0}"/>
      </ParameterList>
    </ParameterList>     





Coordinator
############
 Simulation controller and top-level driver

In the `"cycle driver`" sublist, the user specifies global control of the
simulation, including starting and ending times and restart options.

.. _coordinator-spec:
.. admonition:: coordinator-spec

    * `"start time`" ``[double]`` **0.** Specifies the start of time in model time.
    * `"start time units`" ``[string]`` **"s"** One of "s", "d", or "yr"

    ONE OF

    * `"end time`" ``[double]`` Specifies the end of the simulation in model time.
    * `"end time units`" ``[string]`` **"s"** One of `"s`", `"d`", or `"yr`"

    OR

    * `"end cycle`" ``[int]`` **optional** If provided, specifies the end of the
      simulation in timestep cycles.

      END

    * `"subcycled timestep`" ``[bool]`` **false**  If true, this coordinator creates
      a third State object to store intermediate solutions, allowing for failed
      steps.
    * `"restart from checkpoint file`" ``[string]`` **optional** If provided,
      specifies a path to the checkpoint file to continue a stopped simulation.
    * `"wallclock duration [hrs]`" ``[double]`` **optional** After this time, the
      simulation will checkpoint and end.
    * `"required times`" ``[io-event-spec]`` **optional** An IOEvent_ spec that
      sets a collection of times/cycles at which the simulation is guaranteed to
      hit exactly.  This is useful for situations such as where data is provided at
      a regular interval, and interpolation error related to that data is to be
      minimized.
    * `"PK tree`" ``[pk-typed-spec-list]`` List of length one, the top level
      PK_ spec.

Note: Either `"end cycle`" or `"end time`" are required, and if
both are present, the simulation will stop with whichever arrives
first.  An `"end cycle`" is commonly used to ensure that, in the case
of a time step crash, we do not continue on forever spewing output.

Example:

.. code-block:: xml

   <ParameterList name="cycle driver">
     <Parameter  name="end cycle" type="int" value="6000"/>
     <Parameter  name="start time" type="double" value="0."/>
     <Parameter  name="start time units" type="string" value="s"/>
     <Parameter  name="end time" type="double" value="1"/>
     <Parameter  name="end time units" type="string" value="yr"/>
     <ParameterList name="required times">
       <Parameter name="start period stop" type="Array(double)" value="{0,-1,86400}" />
     </ParameterList>
     <ParameterList name="PK tree">
       <ParameterList name="my richards pk">
         <Parameter name="PK type" type="string" value="richards" />
       </ParameterList>
     </ParameterList>
   </ParameterList>



   

Visualization
##############
 Manages simulation output to disk.

A user may request periodic writes of field data for the purposes of
visualization in the `"visualization`" sublists.

ATS accepts a visualization list for each domain/mesh, including surface and
column meshes.  These are in separate ParameterLists, entitled
`"visualization`" for the main mesh, and `"visualization surface`" on the
surface mesh.  It is expected that, for any addition meshes, each will have a
domain name and therefore admit a spec of the form: `"visualization
DOMAIN-NAME`".

.. _visualization-spec:
.. admonition:: visualization-spec

    * `"file name base`" ``[string]`` **visdump_DOMAIN_data**
    * `"dynamic mesh`" ``[bool]`` **false** Write mesh data for every
      visualization dump; this facilitates visualizing deforming meshes.
    * `"time unit`" ``[string]`` **s** A valid time unit to convert time
      into for output files.  One of `"s`", `"d`", `"y`", or `"yr 365`"

    INCLUDES:
    - ``[io-event-spec]`` An IOEvent_ spec


Example:

.. code-block:: xml

  <ParameterList name="visualization">
    <Parameter name="file name base" type="string" value="visdump_data"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{{0, 100, -1}}" />
    <Parameter name="cycles" type="Array(int)" value="{{999, 1001}}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{{0.0, 10.0, 100.0}}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{{100.0, 25.0, -1.0}}"/>
    <Parameter name="times" type="Array(double)" value="{{101.0, 303.0, 422.0}}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>
  </ParameterList>





Checkpoint
##############
 Manages checkpoint/restart capability.

A user may request periodic dumps of ATS Checkpoint Data in the
`"checkpoint`" sublist.  The user has no explicit control over the
content of these files, but has the guarantee that the ATS run will be
reproducible (with accuracies determined by machine round errors and
randomness due to execution in a parallel computing environment).
Therefore, output controls for Checkpoint Data are limited to file
name generation and writing frequency, by numerical cycle number.
Unlike `"visualization`", there is only one `"checkpoint`" list for
all domains/meshes.

.. _checkpoint-spec:
.. admonition:: checkpoint-spec

    * `"file name base`" ``[string]`` **"checkpoint"**
    * `"file name digits`" ``[int]`` **5**
    * `"single file checkpoint`" ``[bool]`` **true** If true, writes all
      checkpoint to one file.  If false, uses a subdirectory with one file per
      mesh.  false is required if meshes exist on other communicators than
      MPI_COMM_WORLD, but this is toggled if the code detects that this is
      necessary.

    INCLUDES:
    - ``[io-event-spec]`` An IOEvent_ spec

Example:

.. code-block:: xml

  <ParameterList name="checkpoint">
    <Parameter name="cycles start period stop" type="Array(int)" value="{{0, 100, -1}}" />
    <Parameter name="cycles" type="Array(int)" value="{{999, 1001}}" />
    <Parameter name="times start period stop 0" type="Array(double)" value="{{0.0, 10.0, 100.0}}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{{100.0, 25.0, -1.0}}"/>
    <Parameter name="times" type="Array(double)" value="{{101.0, 303.0, 422.0}}"/>
  </ParameterList>

In this example, checkpoint files are written when the cycle number is
a multiple of 100, every 10 seconds for the first 100 seconds, and
every 25 seconds thereafter, along with times 101, 303, and 422.  Files will be written in the form: `"checkpoint00000.h5`".


  


 
Observation
##############
.. _`UnstructuredObservation`:
  Collection of Observations on an unstructured mesh to be written to a common file.


.. _observation-spec:
.. admonition:: observation-spec

    * `"observation output filename`" ``[string]`` user-defined name for the file
      that the observations are written to.
    * `"delimiter`" ``[string]`` **\, ** Delimiter to split columns of the file
    *  `"write interval`" ``[int]`` **1** Interval of observations on which to flush IO files.
    * `"time units`" ``[string]`` **s** Controls the unit of the time column in the observation file.
    * `"domain`" ``[string]`` **""** Can be used to set the communicator which writes (usually not needed).
    * `"observed quantities`" ``[observable-spec-list]`` A list of Observable_
       objects that are all put in the same file.

    INCLUDES:
    - ``[io-event-spec]`` An IOEvent_ spec

Note, for backwards compatibility, an ``observable-spec`` may be directly
included within the `observation-spec` if it is the only variable to be
observed in this file.




.. _`Observable`:
 Collects, reduces, and writes observations during a simulation.

Observations are a localized-in-space but frequent-in-time view of simulation
output, designed to get at useful diagnostic quantities such as hydrographs,
total water content, quantities at a point, etc.  These allow frequent
collection in time without saving huge numbers of visualization files to do
postprocessing.  In fact, these should be though of as orthogonal data queries
to visualization -- vis is pointwise in time but complete in space, while
observations are pointwise/finite in space but complete in time.

Observables describe what is observed -- the variable, region,
integration/averaging scheme, etc.

A user may request any number of specific observables.  Each observable spec
involves a field quantity, a functional reduction operator, a region from which
it will extract its source data, and a list of discrete times for its
evaluation.  The observations are evaluated during the simulation and written
to disk by the UnstructuredObservation_ object.

.. _observable-spec:
.. admonition:: observable-spec

    * `"variable`" ``[string]`` any ATS variable used by any PK, e.g. `"pressure`"
      or `"surface-water_content`"

    * `"region`" ``[string]`` the label of a user-defined region

    * `"location name`" ``[string]`` **cell** The mesh location of the thing to
      be measured, i.e. `"cell`", `"face`", or `"node`"

    * `"number of vectors`" ``[int]`` **1** Number of vector components to write.

    * `"degree of freedom`" ``[int]`` **-1** Degree of freedom to write.  Default
        of -1 implies writing all degrees of freedom.

    * `"functional`" ``[string]`` the type of function to apply to the variable
      on the region.  One of:

      - `"point`" returns the value of the field quantity at a
        point.  The region and location name should result in a single entity being
        selected.

      - `"average`" returns the volume-weighted average of
        the field across all entities in the region selected.  This is likely
        what you want for intensive state variables such as `"temperature`" or
        `"saturation_liquid`".

      - `"integral`" returns the volume-weighted sum of a
        variable over the region.  This should be used for example on intensive
        sources, for instance `"surface-precipitation`", to get the total
        source/sink.

      - `"extensive integral`" returns the sum of an variable
        over the region.  This should be used for extensive quantities such as
        `"water_content`" or `"energy`" which already include the volume in
        their value.

      - `"minimum`" returns the min value over the region

      - `"maximum`" returns the max value over the region

    * `"direction normalized flux`" ``[bool]`` **false** For flux observations,
      dots the face-normal flux with a vector to ensure fluxes are integrated
      pointing the same direction.

    * `"direction normalized flux direction`" ``[Array(double)]`` **optional** For
      flux observations, provides the vector to dot the face normal with.  If this
      is not provided, then it is assumed that the faces integrated over are all
      boundary faces and that the default vector is the outward normal direction
      for each face.

    * `"direction normalized flux relative to region`" ``[string]`` **optional**
      If provided, the flux observation is assumed to be on a set of faces
      which are the exterior of a volumetric region.  This region provides that
      volume, and fluxes are oriented in the "outward normal" direction
      relative to this region's interior.

    * `"time integrated`" ``[bool]`` **false** If true, observe the
      time-integral, observing on all cycles and accumulating the
      backwards-Euler product of dt times the observable at the new time.





PK
###
 The interface for a Process Kernel, an equation or system of equations.

A process kernel represents a single or system of partial/ordinary
differential equation(s) or conservation law(s), and is used as the
fundamental unit for coupling strategies.

Implementations of this interface typically are either an MPC
(multi-process coupler) whose job is to heirarchically couple several
other PKs and represent the system of equations, or a Physical PK,
which represents a single equation.

Note there are two PK specs -- the first is the "typed" spec, which appears in
the "cycle driver" list in the PK tree.  The second is the spec for the base
class PK, which is inherited and included by each actual PK, and lives in the
"PKs" sublist of "main".

.. _pk-typed-spec:
.. admonition:: pk-typed-spec

    * `"PK type`" ``[string]`` One of the registered PK types
    * `"verbose object`" ``[verbose-object-spec]`` **optional** See `Verbose Object`_

.. _pk-spec:
.. admonition:: pk-spec

    * `"PK type`" ``[string]`` One of the registered PK types.  Note this must
      match the corresponding entry in the ``[pk-typed-spec]``
    * `"verbose object`" ``[verbose-object-spec]`` **optional** See `Verbose Object`_

Example:

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="my cool PK">
      <Parameter name="PK type" type="string" value="my cool PK"/>
       ...
    </ParameterList>
  </ParameterList>

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="Top level MPC">
      <Parameter name="PK type" type="string" value="strong MPC"/>
      <ParameterList name="sub PKs">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>




.. contents:: **List of PKs**
   :local:

Base PKs
========
There are several types of PKs, and each PK has its own valid input
spec.  However, there are three main types of PKs, from which nearly
all PKs derive.  Note that none of these are true PKs and cannot stand
alone.

PK: Physical
------------
 A base class with default implementations of methods for a leaf of the PK tree (a conservation equation, or similar).

`PKPhysicalBase` is a base class providing some functionality for PKs which
are defined on a single mesh, and represent a single process model.  Typically
all leaves of the PK tree will inherit from `PKPhysicalBase`.

.. _pk-physical-default-spec:
.. admonition:: pk-physical-default-spec

    * `"domain name`" ``[string]`` Name from the Mesh_ list on which this PK is defined.

    * `"primary variable key`" ``[string]`` The primary variable
      e.g. `"pressure`", or `"temperature`".  Most PKs supply sane defaults.

    * `"initial condition`" ``[initial-conditions-spec]``  See InitialConditions_.

    * `"max valid change`" ``[double]`` **-1** Sets a limiter on what is a
      valid change in a single timestep.  Changes larger than this are declared
      invalid and the timestep shrinks.  By default, any change is valid.
      Units are the same as the primary variable.

    INCLUDES:

    - ``[pk-spec]`` This *is a* PK_.
    - ``[debugger-spec]`` Uses a Debugger_





PK: BDF
-------
 A base class with default implementations of methods for a PK that can be implicitly integrated in time.

`PKBDFBase` is a base class from which PKs that want to use the `BDF`
series of implicit time integrators must derive.  It specifies both the
`BDFFnBase` interface and implements some basic functionality for `BDF`
PKs.

.. _pk-bdf-default-spec:
.. admonition:: pk-bdf-default-spec

    * `"initial time step`" ``[double]`` **1.** Initial time step size [s]

    * `"assemble preconditioner`" ``[bool]`` **true** A flag, typically not set
      by user but by an MPC.

    * `"time integrator`" ``[bdf1-ti-spec]`` **optional**
      A TimeIntegrator_.  Note that this is only required if this PK is not
      strongly coupled to other PKs.

    * `"inverse`" ``[inverse-typed-spec]`` **optional** A Preconditioner_.
      Note that this is only used if this PK is not strongly coupled to other PKs.

    INCLUDES:

    - ``[pk-spec]`` This *is a* PK_.





PK: Physical and BDF
--------------------
 Default implementation of both BDF and Physical PKs.

A base class for all PKs that are both physical, in the sense that they
implement an equation and are not couplers, and BDF, in the sense that they
support the implicit integration interface.  This largely just supplies a
default error norm based on error in conservation relative to the extent of the
conserved quantity.

By default, the error norm used by solvers is given by:

:math:`ENORM(u, du) = |du| / ( a_tol + r_tol * |u| )`

The defaults here are typically good, or else good defaults are set in the
code, so usually are not supplied by the user.


.. _pk-physical-bdf-default-spec:
.. admonition:: pk-physical-bdf-default-spec

    * `"conserved quantity key`" ``[string]`` Name of the conserved quantity.
      Usually a sane default is set by the PK.

    * `"absolute error tolerance`" ``[double]`` **1.0** Absolute tolerance,
      :math:`a_tol` in the equation above.  Unit are the same as the conserved
      quantity.  Note that this default is often overridden by PKs with more
      physical values, and very rarely are these set by the user.

    * `"relative error tolerance`" ``[double]`` **1.0** Relative tolerance,
      :math:`r_tol` in the equation above.  ``[-]`` Note that this default can
      be overridden by PKs with more physical values, and very rarely are these
      set by the user.

    * `"flux error tolerance`" ``[double]`` **1.0** Relative tolerance on the
      flux.  Note that this default is often overridden by PKs with more physical
      values, and very rarely are these set by the user.

    INCLUDES:

    - ``[pk-bdf-default-spec]`` *Is a* `PK: BDF`_
    - ``[pk-physical-default-spec]`` *Is a* `PK: Physical`_





Physical PKs
============
Physical PKs are the physical capability implemented within ATS.

Flow PKs
--------

Flow PKs describe the conservation of mass of water as it flows both
above and below-ground.  Subsurface flow PKs are based on 3D Richards
equation, which describes variably saturated flow in porous media.
Minor variations to this include the incorporation of freeze-thaw
processes.  Surface flow PKs are based on a diffusion wave equation
and Manning's model for sheet flow.  Variations to this also include
the incorporation of freeze-thaw processes.  Finally we include in
flow a "snow distribution" algorithm which takes as input
precipitation and applies it based on the existing surface level
(elevation + water + snowpack), thereby "filling in" low-lying areas
preferentially.  This makes for more accurate snowpacks at fine
scales.

Richards PK
^^^^^^^^^^^
 Two-phase, variable density Richards equation.

Solves Richards equation:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla \cdot \frac{k_r n_l}{\mu} K ( \nabla p + \rho g \hat{z} ) = Q_w

.. _richards-spec:
.. admonition:: richards-spec

   * `"domain`" ``[string]`` **"domain"**  Defaults to the subsurface mesh.

   * `"primary variable key`" ``[string]`` The primary variable associated with
      this PK, typically `"DOMAIN-pressure`" Note there is no default -- this
      must be provided by the user.

   * `"boundary conditions`" ``[list]`` Defaults to Neuman,
      0 normal flux.  See `Flow-specific Boundary Conditions`_

   * `"permeability type`" ``[string]`` **scalar** This controls the number of
      values needed to specify the absolute permeability.  One of:

      - `"scalar`" Requires one scalar value.
      - `"horizontal and vertical`" Requires two values, horizontal then vertical.
      - `"diagonal tensor`" Requires dim values: {xx, yy} or {xx, yy, zz}
      - `"full tensor`". (Note symmetry is required.)  Either {xx, yy, xy} or {xx,yy,zz,xy,xz,yz}.

   * `"water retention evaluator`" ``[wrm-evaluator-spec]`` The water retention
      curve.  This needs to go away, and should get moved to State.

   IF
   * `"source term`" ``[bool]`` **false** Is there a source term?

   THEN
   * `"source key`" ``[string]`` **DOMAIN-water_source** Typically
      not set, as the default is good. ``[mol s^-1]``
   * `"source term is differentiable`" ``[bool]`` **true** Can the source term
      be differentiated with respect to the primary variable?
   * `"explicit source term`" ``[bool]`` **false** Apply the source term from
      the previous time step.
   END

   Math and solver algorithm options:

   * `"diffusion`" ``[pde-diffusion-spec]`` The (forward) diffusion operator,
      see PDE_Diffusion_.

   * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` **optional** The
      inverse of the diffusion operator.  See PDE_Diffusion_.  Typically this
      is only needed to set Jacobian options, as all others probably should
      match those in `"diffusion`", and default to those values.

   * `"surface rel perm strategy`" ``[string]`` **none** Approach for
      specifying the relative permeabiilty on the surface face.  `"clobber`" is
      frequently used for cases where a surface rel perm will be provided.  One
      of:

      - `"none`" : use the upwind direction to determine whether to use the
        boundary face or internal cell
      - `"clobber`" : always use the boundary face rel perm
      - `"max`" : use the max of the boundary face and internal cell values
      - `"unsaturated`" : Uses the boundary face when the internal cell is not
        saturated.

   * `"relative permeability method`" ``[string]`` **upwind with Darcy flux**
      Relative permeability is defined on cells, but must be calculated on
      faces to multiply a flux.  There are several methods commonly used.  Note
      these can significantly change answers -- you don't want to change these
      unless you know what they mena.  One of:

      - `"upwind with Darcy flux`" First-order upwind method that is most common
      - `"upwind with gravity`" Upwinds according to the gravitational flux direction
      - `"cell centered`" This corresponds to the harmonic mean, and is most
        accurate if the problem is always wet, but has issues when it is dry.
      - `"arithmetic mean`" Face value is the mean of the neighboring cells.
        Not a good method.

   Globalization and other process-based hacks:

   * `"modify predictor with consistent faces`" ``[bool]`` **false** In a
      face+cell diffusion discretization, this modifies the predictor to make
      sure that faces, which are a DAE, are consistent with the predicted cells
      (i.e. face fluxes from each sides match).

   * `"modify predictor for flux BCs`" ``[bool]`` **false** Infiltration into
      dry ground can be hard on solvers -- this tries to do the local nonlinear
      problem to ensure that face pressures are consistent with the
      prescribed flux in a predictor.

   * `"modify predictor via water content`" ``[bool]`` **false** Modifies the
      predictor using the method of Krabbenhoft [??] paper.  Effectively does a
      change of variables, extrapolating not in pressure but in water content,
      then takes the smaller of the two extrapolants.

   * `"max valid change in saturation in a time step [-]`" ``[double]`` **-1**
      Rejects timesteps whose max saturation change is greater than this value.
      This can be useful to ensure temporally resolved solutions.  Usually a
      good value is 0.1 or 0.2.

   * `"max valid change in ice saturation in a time step [-]`" ``[double]``
      **-1** Rejects timesteps whose max ice saturation change is greater than
      this value.  This can be useful to ensure temporally resolved solutions.
      Usually a good value is 0.1 or 0.2.

   * `"limit correction to pressure change [Pa]`" ``[double]`` **-1** If > 0,
      this limits an iterate's max pressure change to this value.  Not usually
      helpful.

   * `"limit correction to pressure change when crossing atmospheric [Pa]`" ``[double]``
      **-1** If > 0, this limits an iterate's max pressure change
      to this value when they cross atmospheric pressure.  Not usually helpful.

   Discretization / operators / solver controls:

   * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` **optional**
      The inverse of the accumulation operator.  See PDE_Accumulation_.
      Typically not provided by users, as defaults are correct.

   * `"absolute error tolerance`" ``[double]`` **2750.0** ``[mol]``

   * `"compute boundary values`" ``[bool]`` **false** Used to include boundary
      face unknowns on discretizations that are cell-only (e.g. FV).  This can
      be useful for surface flow or other wierd boundary conditions.  Usually
      provided by MPCs that need them.

   Physics control:

   * `"permeability rescaling`" ``[double]`` **1e7** Typically 1e7 or order
      :math:`sqrt(K)` is about right.  This rescales things to stop from
      multiplying by small numbers (permeability) and then by large number
      (:math:`\rho / \mu`).

   IF
   * `"coupled to surface via flux`" ``[bool]`` **false** If true, apply
      surface boundary conditions from an exchange flux.  Note, if this is a
      coupled problem, it is probably set by the MPC.  No need for a user to
      set it.

   THEN
   * `"surface-subsurface flux key`" ``[string]`` **DOMAIN-surface_subsurface_flux**
   END

   * `"coupled to surface via head`" ``[bool]`` **false** If true, apply
      surface boundary conditions from the surface pressure (Dirichlet).




Permafrost Flow PK
^^^^^^^^^^^^^^^^^^
 A three-phase, thermal Richard's equation with water, water vapor, and ice for permafrost applications.

Note that the only difference between permafrost and richards is in
constitutive relations -- the WRM changes to provide three saturations,
while the water content changes to account for water in ice phase.  As these
are now drop-in field evaluators, there is very little to change in the PK.

In the future, this should not even need a different PK.

.. _permafrost-spec:
.. admonition:: permafrost-spec
  
    * `"saturation ice key`" ``[string]`` **"DOMAIN-saturation_ice"** volume fraction of the ice phase (only when relevant) ``[-]`` Typically the default is correct.

    INCLUDES:

    - ``[richards-spec]`` See `Richards PK`_



Overland Flow PK
^^^^^^^^^^^^^^^^
 Overland flow using the diffusion wave equation.

Solves the diffusion wave equation for overland flow with pressure as a primary variable:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla n_l k \nabla h(p) = Q_w


.. _overland-pressure-spec:
.. admonition:: overland-pressure-spec

    Keys name variables:

    * `"domain`" ``[string]`` **"surface"**  Defaults to the extracted surface mesh.

    * `"primary variable`" ``[string]`` The primary variable associated with
      this PK, typically `"DOMAIN-pressure`" Note there is no default -- this
      must be provided by the user.

    * `"boundary conditions`" ``[list]`` Defaults to Neuman, 0 normal flux.

    * `"overland conductivity evaluator`" ``[list]``
      See `Overland Conductivity Evaluator`_.

    IF

    * `"source term`" ``[bool]`` **false** Is there a source term?

    THEN

    * `"source key`" ``[string]`` **DOMAIN-water_source** Typically
      not set, as the default is good. ``[m s^-1]`` or ``[mol s^-1]``
    * `"water source in meters`" ``[bool]`` **true** Is the source term in ``[m s^-1]``?
    * `"source term is differentiable`" ``[bool]`` **true** Can the source term
      be differentiated with respect to the primary variable?
    * `"explicit source term`" ``[bool]`` **false** Apply the source term from
      the previous time step.

    END

    Math and solver algorithm options:

    * `"diffusion`" ``[pde-diffusion-spec]`` The (forward) diffusion operator,
      see PDE_Diffusion_.

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` **optional** The
      inverse of the diffusion operator.  See PDE_Diffusion_.  Typically this
      is only needed to set Jacobian options, as all others probably should
      match those in `"diffusion`", and default to those values.

    * `"absolute error tolerance`" ``[double]`` **550.** Defaults to 1 cm of
      water.  A small, but significant, amount of water.

    * `"limit correction to pressure change [Pa]`" ``[double]`` **-1** If > 0,
      this limits an iterate's max pressure change to this value.  Not usually
      helpful.

    * `"limit correction to pressure change when crossing atmospheric [Pa]`" ``[double]``
      **-1** If > 0, this limits an iterate's max pressure change
      to this value when they cross atmospheric pressure.  Not usually helpful.

    * `"allow no negative ponded depths`" ``[bool]`` **false** Modifies all
      correction updates to ensure only positive ponded depth is allowed.

    * `"min ponded depth for velocity calculation`" ``[double]`` **1.e-2** For
      ponded depth below this height, declare the velocity 0.

    * `"min ponded depth for tidal bc`" ``[double]`` **0.02** Control on the
      tidal boundary condition.  TODO: This should live in the BC spec?

    INCLUDES:

    - ``[pk-physical-bdf-default-spec]`` A `PK: Physical and BDF`_ spec.

    Everything below this point is usually not provided by the user, but are
    documented here for completeness.

    Keys name variables:

    * `"conserved quantity key`" ``[string]`` **DOMAIN-water_content** Typically
      not set, as the default is good. ``[mol]``
    * `"elevation key`" ``[string]`` **DOMAIN-elevation** Typically
      not set, as the default is good. ``[mol]``
    * `"slope magnitude key`" ``[string]`` **DOMAIN-slope_magnitude** Typically
      not set, as the default is good. ``[mol]``

    Algorithmic parameters:

    * `"coupled to subsurface via flux`" ``[bool]`` **false** Set by MPC.
    * `"coupled to subsurface via head`" ``[bool]`` **false** Set by MPC.

    * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` **optional**
      The inverse of the accumulation operator.  See PDE_Accumulation_.
      Typically not provided by users, as defaults are correct.

    EVALUATORS:

    - `"conserved quantity`"
    - `"water content`"
    - `"cell volume`"
    - `"surface_subsurface_flux`"
    - `"elevation`"
    - `"slope magnitude`"
    - `"overland_conductivity`"
    - `"ponded_depth`"
    - `"pres_elev`"
    - `"source`"


.. todo:
    Nearly all variable name roots are hard-coded here, this should get updated.




Overland Flow with Ice
^^^^^^^^^^^^^^^^^^^^^^
 Two-phase overland flow equation.

This modifies the diffusion wave equation for overland flow that includes
freeze-thaw processes.  This class could completely go away, but it does some
error checking on the input file to make sure freeze-thaw processes are done
correctly.  In the future this service should be done by a preprocessor
generating the input file, and this class would go away.

.. _icy-overland-spec:
.. admonition:: icy-overland-spec

    INCLUDES:

    - ``[overland-pressure-spec]`` See `Overland Flow PK`_.




Snow Distribution PK
^^^^^^^^^^^^^^^^^^^^
 Preferential distribution of snow precip in low-lying areas.

This PK is a heuristic PK that distributes incoming snow precipitation using a
diffusion wave equation.  Think of it as an analogue to overland flow -- it
effectively ensures that new snow "flows downhill," due to a uniformly random
direction and strength wind, and lands on the lowest lying areas.

Tweaking the snow-manning_coefficient lets you play with how uniform the snow
layer ends up.  Most of the parameters are set by your snow precipitation input
data interval.  The details of this are a bit tricky mathematically, and it may
take some fiddling with parameters to do this correctly if your data is not
daily (which all defaults are set for).

.. _snow-distribution-spec:
.. admonition:: snow-distribution-spec

    * `"distribution time`" ``[double]`` **86400.** Interval of snow precip input dataset. `[s]`
    * `"precipitation function`" ``[function-spec]`` Snow precipitation Function_ spec.

    * `"diffusion`" ``[pde-diffusion-spec]`` Diffusion drives the distribution.
      Typically we use finite volume here.  See PDE_Diffusion_

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` Inverse of the
      above.  Likely only Jacobian term options are needed here, as the others
      default to the same as the `"diffusion`" list.  See PDE_Diffusion_.

    * `"inverse`" ``[inverse-typed-spec]`` Inverse_ method for the solve.
    
    Not typically provided by the user, defaults are good:

    * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` See PDE_Accumulation_.
    

.. todo::
    For this PK, all variable root names are hard-coded.  This should get changed.
    




Transport PK
------------

The Transport PK describes the conservation of mass of components transported
with water as it flows. The transport PK is based on the advection-diffusion 
equation, applies to one or more components that are dissolved in the aqueous 
phase, and is currently used in both surface and subsurface compartments. 
The key difference between surface and subsurface transport is in capturing 
the volume of water. In the subsurface, the volume of water is set by the 
porosity and saturation of the porous medium, while in the surface it is set 
by the ponded depth.



Energy PKs
-----------

Energy PKs describe the conservation of energy as it is advected and
diffuses both above and below-ground.  Both surface and subsurface
energy equations are based on a simple advection-diffusion equation,
and include variants with and without freeze-thaw processes.

Energy Base PK
^^^^^^^^^^^^^^
 An advection-diffusion equation for energy.

Solves an advection-diffusion equation for energy:

.. math::
    \frac{\partial E}{\partial t} - \nabla \cdot \kappa \nabla T + \nabla \cdot \mathbf{q} e(T) = Q_w e(T) + Q_e

.. todo:: Document the energy error norm!

.. _energy-pk-spec:
.. admonition:: energy-pk-spec

    * `"domain`" ``[string]`` **"domain"**  Defaults to the subsurface mesh.

    * `"primary variable`" ``[string]`` The primary variable associated with
      this PK, typically `"DOMAIN-temperature`" Note there is no default -- this
      must be provided by the user.

    * `"boundary conditions`" ``[list]`` Defaults to 0 diffusive flux
      boundary condition.  See `Energy-specific Boundary Conditions`_

    * `"thermal conductivity evaluator`" ``[list]``
      The thermal conductivity.  This
      needs to go away, and should get moved to State.

    * `"absolute error tolerance`" ``[double]`` **76.e-6** A small amount of
      energy, see error norm. `[MJ]`

    * `"upwind conductivity method`" ``[string]`` **arithmetic mean** Method of
      moving cell-based thermal conductivities onto faces.  One of:

      - `"arithmetic mean`" the default, average of neighboring cells
      - `"cell centered`" harmonic mean

    IF

    * `"explicit advection`" ``[bool]`` **false** Treat the advection term implicitly.

    ELSE

    * `"supress advective terms in preconditioner`" ``[bool]`` **false**
      Typically subsurface energy equations are strongly diffusion dominated,
      and the advective terms may add little.  With this flag on, we ignore
      theem in the preconditioner, making an easier linear solve and often not
      negatively impacting the nonlinear solve.

    * `"advection preconditioner`" ``[list]`` **optional**
      Typically defaults are correct.

    END

    * `"diffusion`" ``[pde-diffusion-spec]`` See PDE_Diffusion_, the diffusion operator.

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` See
      PDE_Diffusion_, the inverse operator.  Typically only adds Jacobian
      terms, as all the rest default to those values from `"diffusion`".

    IF

    * `"source term`" ``[bool]`` **false** Is there a source term?

    THEN

    * `"source key`" ``[string]`` **DOMAIN-total_energy_source** Typically
      not set, as the default is good. ``[MJ s^-1]``

    * `"source term is differentiable`" ``[bool]`` **true** Can the source term
      be differentiated with respect to the primary variable?

    * `"source term finite difference`" ``[bool]`` **false** If the source term
      is not diffferentiable, we can do a finite difference approximation of
      this derivative anyway.  This is useful for difficult-to-differentiate
      terms like a surface energy balance, which includes many terms.

    END

    Globalization:

    * `"modify predictor with consistent faces`" ``[bool]`` **false** In a
      face+cell diffusion discretization, this modifies the predictor to make
      sure that faces, which are a DAE, are consistent with the predicted cells
      (i.e. face fluxes from each sides match).

    * `"modify predictor for freezing`" ``[bool]`` **false** A simple limiter
      that keeps temperature corrections from jumping over the phase change.

    * `"limit correction to temperature change [K]`" ``[double]`` **-1.0** If >
      0, stops nonlinear updates from being too big through clipping.

    The following are rarely set by the user, as the defaults are typically right.

    * `"advection`" ``[list]`` **optional** The PDE_Advection_ spec.  Only one
      current implementation, so defaults are typically fine.

    * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` **optional**
      The inverse of the accumulation operator.  See PDE_Accumulation_.
      Typically not provided by users, as defaults are correct.

    IF
    * `"coupled to surface via flux`" ``[bool]`` **false** If true, apply
      surface boundary conditions from an exchange flux.  Note, if this is a
      coupled problem, it is probably set by the MPC.  No need for a user to
      set it.
    THEN
    * `"surface-subsurface energy flux key`" ``[string]`` **DOMAIN-surface_subsurface_energy_flux**
    END

    * `"coupled to surface via temperature`" ``[bool]`` **false** If true, apply
      surface boundary conditions from the surface temperature (Dirichlet).

    KEYS:
    - `"conserved quantity`" **DOMAIN-energy** The total energy :math:`E` `[MJ]`
    - `"energy`" **DOMAIN-energy** The total energy :math:`E`, also the conserved quantity. `[MJ]`
    - `"water content`" **DOMAIN-water_content** The total mass :math:`\Theta`, used in error norm `[mol]`
    - `"enthalpy`" **DOMAIN-enthalpy** The specific enthalpy :math`e` `[MJ mol^-1]`
    - `"flux`" **DOMAIN-mass_flux** The mass flux :math:`\mathbf{q}` used in advection. `[mol s^-1]`
    - `"diffusive energy`" **DOMAIN-diffusive_energy_flux** :math:`\mathbf{q_e}` `[MJ s^-1]`
    - `"advected energy`" **DOMAIN-advected_energy_flux** :math:`\mathbf{q_e^{adv}} = q e` `[MJ s^-1]`
    - `"thermal conductivity`" **DOMAIN-thermal_conductivity** Thermal conductivity on cells `[W m^-1 K^-1]`
    - `"upwinded thermal conductivity`" **DOMAIN-upwinded_thermal_conductivity** Thermal conductivity on faces `[W m^-1 K^-1]`

    EVALUATORS:
    - `"source term`" **optional** If source key is provided.
    - `"enthalpy`"
    - `"cell volume`"
    - `"thermal conductivity`"
    - `"conserved quantity`"
    - `"energy`"




Two-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 An advection-diffusion equation for energy in two phases.

This is simply a subsurface energy equation that places a few more requirements
on the base class.  It could probably go away if we refactor to remove
hard-coded evaluators.

.. _energy-two-phase-pk-spec:
.. admonition:: energy-two-phase-pk-spec

    INCLUDES:

    - ``[energy-pk-spec]``  See `Energy Base PK`_




Three-Phase subsurface Energy PK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 An advection-diffusion equation for energy in three phases.

This is simply a subsurface energy equation that places a few more requirements
on the base class.  It could probably go away if we refactor to remove
hard-coded evaluators.

.. _energy-three-phase-pk-spec:
.. admonition:: energy-three-phase-pk-spec

    INCLUDES:

    - ``[energy-two-phase-pk-spec]`` See  `Two-Phase subsurface Energy PK`_




Overland energy with Ice
^^^^^^^^^^^^^^^^^^^^^^^^
 An advection-diffusion equation for surface energy in two phases.

This is simply a surface energy equation that places a few more requirements
on the base class.  It could probably go away if we refactor to remove
hard-coded evaluators.

.. _energy-surface-ice-pk-spec:
.. admonition:: energy-surface-ice-pk-spec

    These are typically not set by the user:

    * `"coupled to subsurface via temperature`" ``[bool]`` **false** A coupling
      scheme, provided by MPC.

    * `"coupled to subsurface via flux`" ``[bool]`` **false** A coupling
      scheme, provided by MPC.

    * `"subsurface domain name`" ``[string]`` **optional** If one of the above
      coupling schemes is turned on, we need to know the subsurface mesh.
      Provided by MPC.

    INCLUDES:

    - ``[energy-pk-spec]``  See `Energy Base PK`_






Surface Energy Balance PKs
------------------------------

Integrated hydrology is not much use without significant process
complexity in source terms coming from the ecohydrologic environment.
These include straightforward sources, like precipitation, but also
more complicated ones such as evaporation and transpiration.

These terms are almost always tied up in a surface energy balance --
evaporation and transpiration are driven by vapor pressure gradients
between the atmosphere and the surface (either snow, ponded water,
soil, or leaf).  Solving a surface energy balance often requires
providing a bunch of terms, including radiated energy, conducted
energy, latent and sensible heat models, etc.

ATS currently has several approaches to calculating these -- see
`ats-demos <https://github.com/amanzi/ats-demos>`_ examples on
ecohydrology for a more in-depth discussion.

Balance Equation
^^^^^^^^^^^^^^^^
 A simple conservation ODE.

This is a very simple vector of ODEs, useful in balance equations, where the
time derivative of a conserved quantity is determined by a bunch of sources and
sinks.

.. math::
    \frac{\partial \Phi }{\partial t} = \sum_i Q_i

.. _balance-pk-spec:
.. admonition:: balance-pk-spec

    * `"primary variable key`" ``[string]`` The primary variable associated with
      this PK.  Note there is no default -- this must be provided by the user.

    * `"conserved quantity key`" ``[string]`` The conserved quantity :math:`\Phi`

    * `"source key`" ``[string]`` **DOMAIN-source_sink** Units are in conserved
      quantity per second per cell volume.

    * `"time discretization theta`" ``[double]`` **1.0** :math:`\theta` in a
      Crank-Nicholson time integration scheme.  1.0 implies fully implicit, 0.0
      implies explicit, 0.5 implies C-N.

    * `"modify predictor positivity preserving`" ``[bool]`` **false** If true,
      predictors are modified to ensure that the conserved quantity is always > 0.

    * `"absolute error tolerance`" ``[double]`` **550.0** a_tol in the standard
      error norm calculation.  Defaults to a small amount of water.  Units are
      the same as the conserved quantity.

    INCLUDES:

    - ``[pk-physical-bdf-default-spec]``





Snow Balance Equation
^^^^^^^^^^^^^^^^^^^^^
 An implicit PK for surface balance snow SWE conservation.

This is a balance PK whose conserved quantity is snow SWE.  The energy balance
comes in as it provides the energy needed to melt snow.  So source terms
include snow precipitation and snowmelt.  It also manages snow density, which
should get rethought a bit.

There is also some wierd hackiness here about area fractions -- see ATS Issue
#8

.. _subgrid-balance-pk-spec:
.. admonition:: subgrid-balance-pk-spec

    * `"absolute error tolerance`" ``[double]`` **0.01** ``[m]``

    INCLUDES:

    - ``[balance-pk-spec]`` This *is a* `Balance Equation`_

    Not typically set by user, defaults work:

    * `"conserved quantity key`" ``[string]`` **LAYER-snow_water_equivalent**
      Sets the default conserved quantity key, so this is likely not supplied
      by the user. `[m]`
    * `"snow density key`" ``[string]`` **LAYER-density** Default snow density
      key. `[kg m^-3]`
    * `"snow age key`" ``[string]`` **LAYER-age** Default snow age key. `[d]`
    * `"new snow key`" ``[string]`` **LAYER-source** Default new snow key. `[m SWE s^-1]`
    * `"area fractions key`" ``[string]`` **LAYER-fractional_areas** Subgrid
      model fractional areas, see note above. `[-]`
    * `"snow death rate key`" ``[string]`` **LAYER-death_rate** Deals with last
      tiny bit of snowmelt.




Biogeochemistry
---------------

To accurately predict watershed ecohydrology, a carbon cycle model is
needed to predict transpiration.  By simulating a carbon cycle, we are
able to predict the rate of photosynthesis as a function of space and
time, and photosynthesis governs root water uptake.  Currently only
one big-leaf model is available, but ongoing work is wrapping a
generalized Common/Colorado Land Model based on that developed within
the ParFlow team, and another ongoing project is working on wrapping
kernels from E3SM's Land Model.

Biogeochemistry -- Monolithic Version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Above and below-ground carbon cycle model.

This is a multi-leaf layer, big-leaf vegetation model coupled to a Century
model for belowground carbon decomposition.

It leverages a PFT-based structure which allows multiple height-sorted PFTs to
coexist on the same grid cells, with the shorter PFTs getting whatever light is
left in the understory.

The implementation is based on an old, standalone code by Chonggang Xu, and
adapted for ATS.  While this is not simple, it is called BGC simple as it is
about the least amount of complexity required to get a reasonable carbon cycle
into ATS.

Outputs of this include transpiration, a critical sink for hydrology, as it
solves photosynthesis based on water availability.

Note this is an "explicit update PK," or effectively a forward Euler timestep
that is not written in ODE form.

Note this works on both the surface (vegetation) and subsurface (decomposition)
meshes.  **It is required** that the subsurface mesh is a "columnar" mesh, and
that build_columns in the subsurface Mesh_ spec has been supplied.

.. _bgc-simple-spec:
.. admonition:: bgc-simple-spec

  * `"initial time step`" ``[double]`` **1.0** Initial time step size `[s]`

  * `"number of carbon pools`" ``[int]`` **7** Unclear whether this can actually change?

  * `"soil carbon parameters`" ``[soil-carbon-spec-list]`` List of soil carbon parameters by soil mesh partition region name.

  * `"pft parameters`" ``[pft-spec-list]`` List of PFT parameters by PFT name.

  * `"latitude [degrees]`" ``[double]`` **60** Latitude of the simulation in degrees.  Used in radiation balance.

  * `"wind speed reference height [m]`" ``[double]`` **2.0** Reference height of the wind speed dataset.

  * `"cryoturbation mixing coefficient [cm^2/yr]`" ``[double]`` **5.0** Controls diffusion of carbon into the subsurface via cryoturbation.

  * `"leaf biomass initial condition`" ``[initial-conditions-spec]`` Sets the leaf biomass IC.

  * `"domain name`" ``[string]`` **domain**

  * `"surface domain name`" ``[string]`` **surface**

  * `"transpiration key`" ``[string]`` **DOMAIN-transpiration** The distributed transpiration flux `[mol s^-1]`

  * `"shaded shortwave radiation key`" ``[string]``
    **SURFACE_DOMAIN-shaded_shortwave_radiation** Shortwave radiation that gets
    past the canopy and teo the bare ground for soil evaporation. `[W m^-2]`

  * `"total leaf area index key`" ``[string]`` **SURFACE_DOMAIN-total_leaf_area_index** Total LAI across all PFTs.

  EVALUATORS:

  - `"temperature`" The soil temperature `[K]`
  - `"pressure`" soil mafic potential `[Pa]`
  - `"surface-cell_volume`" `[m^2]`
  - `"surface-incoming shortwave radiation`" `[W m^-2]`
  - `"surface-air_temperature`" `[K]`
  - `"surface-relative_humidity`" `[-]`
  - `"surface-wind_speed`" `[m s^-1]`
  - `"surface-co2_concentration`" `[ppm]`
  





Deformation
-------------

The unstructured mesh framework we use provides the opportunity to
include deformation of the mesh.  This deformation can be done in two
ways -- either node coordinate changes are provided, or volumetric
changes are provided, and the code attempts to iterate toward a global
coordinate change that satisfies these volumetric changes.  The latter
can be somewhat fragile for large deformation, but it does allow
simple deformation such as small, somewhat uniform subsidence.  The
volumetric deformation PK below does this based on a volumetric change
given by loss of bulk ice.

Volumetric Deformation
^^^^^^^^^^^^^^^^^^^^^^
 Subsidence through bulk ice loss and cell volumetric change.

This process kernel provides for going from a cell volumetric change to an
updated unstructured mesh, and can be coupled sequentially with flow to solve
problems of flow in a subsiding porous media.

Note that all deformation here is vertical, and we assume that the subsurface
mesh is **perfectly columnar** and that the "build columns" parameter has been
given to the subsurface mesh.  See the Mesh_ spec for more.

The process here is governed through two options, the "deformation mode" and
the "deformation strategy."

The deformation mode describes how the cell volume change is calculated.  There
are three options here:

- "prescribed" uses a function to precribe the volume changes as a function of (t,x,y,z).

- "structural" decreases the cell volume if the porosity is above a prescribed
  "structurally connected matrix" porosity.  Think of this as bulk ice
  "propping up" the soil grains -- as that bulk ice melts, it reduces porosity
  toward the porosity in at which grains start to touch again and can be
  structurally sound.

- "saturation" is a heuristic that considers the liquid saturation directly,
  and tries to relax the liquid saturation back toward a value that is
  consistent with what the thawed soil should be.

.. todo: Move this into an evaluator!

The deformation strategy describes how the cell volume change is turned into
node coordinate changes.  Three options are available:

- "average" simply takes the average of volume change/surface area and
  horizontally averages this quantity across all neighbors.  While this has the
  advantage of being simple, it has issues when thaw gradients in the
  horizontal are not zero, as it may result in the loss of volume in a fully
  frozen cell, blowing up the pressure and breaking the code.  This is great
  when it works, but it almost never works in real problems.

- "global optimization" attempts to directly form and solve the minimization
  problem to find the nodal changes that result in the target volumetric
  changes.  Note this has issues with overfitting, so penalty methods are used
  to smooth the solution of the problem.  This is not particularly robust.

- "mstk implementation" MSTK implements an iterated, local optimization method
  that, one-at-a-time, moves nodes to try and match the volumes.  This has
  fewer issues with overfitting, but doesn't always do sane things, and can be
  expensive if iterations don't work well.  This is not particularly robust
  either, but it seems to be the preferred method for now.

.. todo:: Check that "global optimization" even works? --etc
    
  

.. _volumetric-deformation-pk-spec:
.. admonition:: volumetric-deformation-pk-spec

    * `"max time step [s]`" ``[double]`` **inf** Sets a maximum time step size.

    * `"deformation mode`" ``[string]`` **prescribed** See above for
      descriptions.  One of: `"prescribed`", `"structural`", or `"saturation`".

    * `"deformation strategy`" ``[string]`` **global optimization** See above
      for descriptions.  One of `"average`", `"global optimization`", or `"mstk
      implementation`"

    * `"domain name`" ``[string]`` **domain**  The mesh to deform.

    * `"surface domain name`" ``[string]`` **surface** The surface mesh.

    * `"deformation function`" ``[function-spec]`` **optional** Only used if
      "deformation mode" == "prescribed"

    EVALUATORS:
    - `"saturation_ice`"
    - `"saturation_liquid`"
    - `"saturation_gas`"
    - `"base_porosity`"
    - `"porosity`"
    - `"cell volume`"
      
    INCLUDES:

    - ``[pk-physical-default-spec]``





MPC
===

Multi-process-couplers or MPCs couple other PKs.  They also are PKs
themselves, in that they implement the PK interface.  So MPCs can also
couple other MPCs.  There are a few common "base" MPCs which do the
simplest form of coupling -- sequential and globally implicit (with a
diagonal preconditioner).  Then there are specific couplers which know
more about their coupled sub-PKs, and can do more complicated things
(for instance, adding off-diagonal block entries to the
preconditioner).

MPCs are also used to couple across domains -- for instance integrated
hydrology is a surface+subsurface flow coupler.  They also can do
fancier things like drape a bunch of subgrid columns off of a mesh, or
other things.  Think of these as the custom couplers.

Base MPC
--------
 Multi process coupler base class.

A multi process coupler is a PK (process kernel) which coordinates and couples
several PKs.  Each of these coordinated PKs may be MPCs themselves, or physical
PKs.  Note this does NOT provide a full implementation of PK -- it does not
supply the AdvanceStep() method.  Therefore this class cannot be instantiated, but
must be inherited by derived classes which finish supplying the functionality.
Instead, this provides the data structures and methods (which may be overridden
by derived classes) for managing multiple PKs.

Most of these methods simply loop through the coordinated PKs, calling their
respective methods.

.. _mpc-spec:
.. admonition:: mpc-spec

    * `"PKs order`" ``[Array(string)]`` Provide a specific order to the
      sub-PKs; most methods loop over all sub-PKs, and will call the sub-PK
      method in this order.

    INCLUDES:

    - ``[pk-spec]`` *Is a* PK_.





WeakMPC
-------
 Multi process coupler for sequential coupling.

Noniterative sequential coupling simply calls each PK's AdvanceStep() method in
order.

.. _weak-mpc-spec:
.. admonition:: weak-mpc-spec

    INCLUDES:

    - ``[mpc-spec]`` *Is a* MPC_.




StrongMPC
----------
 Multi process coupler for globally implicit (strong) coupling.

Globally implicit coupling solves all sub-PKs as a single system of equations.  This can be completely automated when all PKs are also `PK: BDF`_ PKs, using a block-diagonal preconditioner where each diagonal block is provided by its own sub-PK.

.. _strong-mpc-spec:
.. admonition:: strong-mpc-spec

    INCLUDES:

    - ``[mpc-spec]`` *Is a* MPC_.
    - ``[pk-bdf-default-spec]`` *Is a* `PK: BDF`_.
    



Physical MPCs
===============

Coupling is an art, and often requires special off-diagonal work for
globally implicit coupling, and fancy games can be played with domains
to couple across domain interfaces both implicitly and sequentially.
Physical MPCs derive from default MPCs to provide special
implementations of some methods.

Coupled Water MPC
-----------------
 A coupler which integrates surface and subsurface flow implicitly.

Couples Richards equation to surface water through continuity of both pressure
and fluxes.  This leverages subsurface discretizations that include face-based
unknowns, and notes that those face unknowns that correspond to surface faces
are co-located with the surface cell pressure, and therefore are equivalent.
In this approach (described in detail in a paper that is in review), the
surface equations are directly assembled into the subsurface discrete operator.

.. _mpc-coupled-water-spec:
.. admonition:: mpc-coupled-water-spec

   * `"PKs order`" ``[Array(string)]`` The use supplies the names of the
     coupled PKs.  The order must be {subsurface_flow_pk, surface_flow_pk}
     (subsurface first).

   * `"subsurface domain name`" ``[string]`` **domain** 

   * `"surface domain name`" ``[string]`` **surface** 

   * `"water delegate`" ``[mpc-delegate-water-spec]`` A `Coupled Water
     Globalization Delegate`_ spec.

   INCLUDES:

   - ``[strong-mpc-spec]`` *Is a* StrongMPC_




Coupled Cells MPC
-----------------
 A coupler which solves two PDEs on the same domain.

This is a StrongMPC which uses a preconditioner in which the
block-diagonal cell-local matrix is dense.  If the system looks something
like:

A( y1, y2, x, t ) = 0
B( y1, y2, x, t ) = 0

where y1,y2 are spatially varying unknowns that are discretized using the MFD
method (and therefore have both cell and face unknowns), an approximation to
the Jacobian is written as

[  dA_c/dy1_c  dA_c/dy1_f   dA_c/dy2_c       0      ]
[  dA_f/dy1_c  dA_f/dy1_f      0              0      ]
[  dB_c/dy1_c     0          dB_c/dy2_c  dB_c/dy2_f ]
[      0           0          dB_f/dy2_c  dB_f/dy2_f ]


Note that the upper left block is the standard preconditioner for the A
system, and the lower right block is the standard precon for the B system,
and we have simply added cell-based couplings, dA_c/dy2_c and dB_c/dy1_c.

Most commonly this is used to couple flow and energy equations on the same
mesh.  In the temperature/pressure system, these extra blocks correspond to

.. math::
    \frac{\partial \Theta}{\partial T} \; , \; \frac{\partial E}{\partial p}

.. _mpc-coupled-cells-spec:
.. admonition:: mpc-coupled-cells-spec

    * `"domain name`" ``[string]`` Domain of simulation
    * `"conserved quantity A`" ``[string]`` Key of the first sub-PK's conserved quantity.
    * `"conserved quantity B`" ``[string]`` Key of the second sub-PK's conserved quantity.
    * `"primary variable A`" ``[string]`` Key of the first sub-PK's primary variable.
    * `"primary variable B`" ``[string]`` Key of the second sub-PK's primary variable.
    * `"no dA/dy2 block`" ``[bool]`` **false** Excludes the dA_c/dy2_c block above.
    * `"no dB/dy1 block`" ``[bool]`` **false** Excludes the dB_c/dy1_c block above.
    
    INCLUDES:

    - ``[strong-mpc-spec]`` *Is a* StrongMPC_.
    



Subsurface MPC
--------------
 A coupler which solves flow and energy in the subsurface.

This MPC provides most nearly all terms for an approximate Jacobian for
coupling three-phase Richards equation (the `Permafrost Flow PK`_) to the
three-phase Energy equation (the `Three-Phase subsurface Energy PK`_).

Many options are provided for turning on and off various aspects of this
Jacobian, so it is useful to mathematically write out these terms.  The
equations are:

.. math::
    \frac{\partial \Theta}{\partial t} - \nabla \frac{k_r n_l}{\mu} K ( \nabla p + \rho g \cdot \hat{z} ) = Q_w \\
    \frac{\partial E}{\partial t} - \nabla \cdot \kappa \nabla T + \nabla \cdot \mathbf{q} e(T) = Q_w e(T) + Q_e

Note that all of the following are dependent on :math:`p` and/or :math:`T`:

.. math::
    \Theta(p,T), k_r(p,T), n_l(p,T), \mu(T), \rho(p,T), E(p,T), \kappa(p,T), e(T)

Also, both source terms :math:`Q_w` and :math:`Q_e` may or may not depend on :math:`p` and :math:`T`.

Note also that the Darcy flux :math:`\mathbf{q}` used in the advection of energy is given by the Darcy flux:

.. math::
    \mathbf{q} = -\frac{k_r n_l}{\mu} K ( \nabla p + \rho g \cdot \rho g \cdot \hat{z} )

Differentiating these two equations in their two unknowns gives the following four blocks in the approximate Jacobian:

:math:`\frac{\partial F_1}{\partial p}`: this is the Richards equation diagonal block, and is controlled inside that PK.

:math:`\frac{\partial F_1}{\partial T}` includes terms for:

- :math:`\frac{\partial \Theta}{\partial T}` This term is the cell-local diagonal block.

- The partial derivative of the divergence of the Darcy flux with respect to
  temperature is dominated by :math:`\frac{\partial}{\partial T} \frac{k_r
  n_l}{\mu}`.  This is because the relative permeability is strongly
  dependent upon phase change (the freezing equals drying approximation).  This
  term is referred to as the "d div q / dT" term.

:math:`\frac{\partial F_2}{\partial p}` includes terms for:

- :math:`\frac{\partial E}{\partial p}` This term is the cell-local diagonal block.

- The partial derivative of the energy diffusion term with respect to pressure
  is dominated by :math:`\frac{\partial \kappa}{\partial p}` through phase
  change -- at a constant temperature, but changing pressure, phase change can
  result in large changes to thermal conductivity.  This is referred to as the
  "div K grad T / dp" term.

:math:`\frac{\partial F_2}{\partial T}`: this is the energy equation diagonal block, and is controlled inside that PK.

Also, at this level, where we know more about the flux used in the energy
equation (it is the Darcy flux), we can do a better approximation of the
derivative of the advection of energy term with respect to both temperature and
pressure.  For instance, enthalpy is only weakly dependent on pressure, so we
can use the derivative of the divergence of the Darcy flux with respect to
pressure (from the Richards block) in the advection term in the
:math:`\frac{\partial F_2}{\partial p}` block, and approximate
:math:`\frac{\partial k_r}{\partial T}` in the advection term as well.  These
terms are referred to as "div hq / dp,T terms".  Note the missing initial "d"
here relative to other terms.

The behavior of this MPC's preconditioner can be set by an option,
`"preconditioner type`".  Really users should not change this from the default,
except in expert cases or for comparison's sake, but the options are:

- `"picard`" is the default, this uses all available terms, and enables the
  "suppress" options for finer-grained control.

- `"none`" No preconditioner never works.

- `"block diagonal`" This is what one would get from the default StrongMPC_.  This probably never works.

- `"no flow coupling`" This keeps the accumulation terms, but turns off all the
  non-local blocks.  This is equivalent to `Coupled Cells MPC`_.

- `"ewc`" **CURRENTLY DEPRECATED/BROKEN/DISABLED** In addition to the
  `"picard`" coupling, this also *always* does a change of variables, whereby
  we first invert to calculate primary variable corrections, then do a change
  of variables to calculate the linearized corrections in energy and water
  content space.  We then apply those corrections, and invert to find the
  primary variable changes that would have made those corrections.  This is
  called the "energy and water content" algorithm, and is related to similar
  variable changing approaches by Krabbenhoft (for flow) and Knoll (for
  energy), but in the multivariate approach.  This is somewhat bad, becuase
  while it fixes some corrections, it breaks others.

- `"smart ewc`" **CURRENTLY DEPRECATED/BROKEN/DISABLED** Does the `"ewc`"
  algorithm above, but tries to be smart about when to do it.  This algorithm
  helps when we are about to fall off of the latent heat cliff.  If we can
  guess when to do it, we have a better chance of not breaking things.  This
  seems like it ought to be helpful, but often doesn't do as much as one might
  hope.


Note this "ewc" algorithm is just as valid, and more useful, in the predictor
(where it is not deprecated/disabled).  There, we extrapolate a change in
pressure and temperature, but often do better to extrapolate in water content
and energy space, then invert (locally) for pressure and temperature
corrections that meet that extrapolation.  Both of these globalization
algorithms are supported by the `EWC Globalization Delegate`_ object.

.. _mpc-subsurface-spec:
.. admonition:: mpc-subsurface-spec

    * `"domain name`" ``[string]`` Domain of simulation

    * `"preconditioner type`" ``[string]`` **picard** See the above for
      detailed descriptions of the choices.  One of: `"none`", `"block
      diagonal`", `"no flow coupling`", `"picard`", `"ewc`", and `"smart ewc`".

    * `"supress Jacobian terms: div hq / dp,T`" ``[bool]`` **false** If using picard or ewc, do not include this block in the preconditioner.
    * `"supress Jacobian terms: d div q / dT`" ``[bool]`` **false** If using picard or ewc, do not include this block in the preconditioner.
    * `"supress Jacobian terms: d div K grad T / dp`" ``[bool]`` **false** If using picard or ewc, do not include this block in the preconditioner.

    * `"ewc delegate`" ``[mpc-delegate-ewc-spec]`` A `EWC Globalization Delegate`_ spec.

    INCLUDES:

    - ``[strong-mpc-spec]`` *Is a* StrongMPC_.

 


Surface MPC
--------------


Permafrost MPC
--------------
 A coupler which solves flow and energy both surface and subsurface.

This MPC handles the coupling of surface energy and flow to subsurface energy
and flow for integrated hydrology with freeze/thaw processes.

.. _mpc-permafrost-spec:
.. admonition:: mpc-permafrost-spec

   * `"PKs order`" ``[Array(string)]`` The user supplies the names of the
     coupled PKs.  The order must be {subsurface_flow_pk, subsurface_energy_pk,
     surface_flow_pk, surface_energy_pk}.

   * `"subsurface domain name`" ``[string]`` **domain** 

   * `"surface domain name`" ``[string]`` **surface** 

   * `"mass exchange flux key`" ``[string]`` **SURFACE_DOMAIN-surface_subsurface_flux**

   * `"energy exchange flux key`" ``[string]`` **SURFACE_DOMAIN-surface_subsurface_energy_flux**

   * `"water delegate`" ``[mpc-delegate-water-spec]`` A `Coupled Water
     Globalization Delegate`_ spec.

   INCLUDES:

   - ``[mpc-subsurface-spec]`` *Is a* `Subsurface MPC`_
    
 



Globalization Delegates
=======================

Globalization is the art of convincing a solver to find the solution.
Remember -- physics typically cares very little about *how* you get to
a solution, only that you get there.  If you can guess or otherwise
find the solution physically, without doing fancy math, go for it!
These delegates are handy utility classes which are used by MPCs to
effeciently leverage physics understanding in the mathematical solvers
to nudge the solver in the direction of a reasonable solution, or to
keep a solver from going off into a part of space which is totally
unphysical.  These can often make the difference between converging
and not converging.

Much of the efficiency of ATS comes from these delegates, and more of
them are always welcome contributions.

Coupled Water Globalization Delegate
------------------------------------
 Globalization for nonlinearity around the appearance/disappearance of surface water.

The water delegate works to deal with discontinuities/strong nonlinearities
when surface cells shift from dry to wet (i.e. the surface pressure goes from <
atmospheric pressure to > atmospheric pressure.

These methods work to alter the predictor around this nonlinearity.

.. _mpc-delegate-water-spec:
.. admonition:: mpc-delegate-water-spec

    * `"modify predictor with heuristic`" ``[bool]`` **false** This simply
      limits the prediction to backtrack to just above atmospheric on both the
      first and second timesteps that take us over atmospheric.

    * `"modify predictor damp and cap the water spurt`" ``[bool]`` **false** The
      second both limits (caps) and damps all surface cells to ensure that all
      nearby cells are also not overshooting.  This is the preferred method.

    These methods work to alter the preconditioned correction for the same
    reasons described above.

    * `"global water face limiter`" ``[double]`` **1.e99** This is simply a limit
      to the maximum allowed size of the correction (in [Pa]) on all faces.  Any
      correction larger than this is set to this.

    * `"cap the water spurt`" ``[bool]`` **false** If a correction takes the
      pressure on a surface cell from below atmospheric (dry) to above (wet),
      the correction is set to a value which results in the new iterate to being
      CAP_SIZE over atmospheric.

    * `"damp the water spurt`" ``[bool]`` **false** A damping factor (less than
      one) is calculated to multiply the correction such that the largest
      correction takes a cell to just above atmospheric.  All faces (globally)
      are affected.

    * `"damp and cap the water spurt`" ``[bool]`` **false** None of the above
      should really be used.  Capping, when the cap is particularly severe,
      results in faces whose values are very out of equilibrium with their
      neighboring cells which are not capped.  Damping results in a tiny
      timestep in which, globally, at MOST one face can go from wet to dry.
      This looks to do a combination, in which all things are damped, but faces
      that are initially expected to go from dry to wet are pre-scaled to
      ensure that, when damped, they are also (like the biggest change) allowed
      to go from dry to wet (so that multiple cells can wet in the same step).
      This is the preferred method.

    In these methods, the following parameters are useful:

    * `"cap over atmospheric`" ``[double]`` **100** This sets the max size over
      atmospheric to which things are capped or damped. `[Pa]`

 


EWC Globalization Delegate
--------------------------
 Globalization for nonlinearity associated with phase change and latent heat.

The EWC delegate works to deal with strong nonlinearities associated with
latent heat and phase change.  Provided a change in primary variables pressure
and temperature, it works by first multiplying those changes by the local
Jacobian matrix, :math:`\frac{\partial \left\{ \Theta, E \right\} }{ \partial
\left\{ p, T \right\} }` to calculate changes in water content and energy, then
calculating the new water content and energy and inverting the functions
:math:`\Theta(p,T), E(p,T)` to determine what pressure and temperature would
have resulted in those values.  This provides a corrected change in the primary
variables.

Conceptually, this is a "more robust" choice in nonlinearities associated with
phase change, where the derivatives go from small to large to small again, and
small changes in pressure and temperature result in large changes in water
content and energy.

This delegate manages these globalization strategies, which can be used both in
modifying the correction supplied by a nonlinear iterate, and in modifying a
predictor, the extrapolated projection (from previous timesteps) that
provides the initial guess to the nonlinear solve.

.. _mpc-delegate-ewc-spec:
.. admonition:: mpc-delegate-ewc-spec

    * `"verbose object`" ``[verbose-object-spec]`` See `Verbose Object`_.
    
    * `"PK name`" ``[string]`` Name of the owning PK -- simply for logging and
      debugging.
    * `"domain name`" ``[string]`` **"domain"** The mesh.

    * `"preconditioner type`" ``[string]`` When to use EWC on the nonlinear
      iterate's correction.  One of:

      - `"none`" Never do EWC
      - `"ewc`" Always do EWC
      - `"smart ewc`" Attempt EWC when it seems likely it will be useful and
        take the EWC correction if it is smaller than the standard correction.

    * `"predictor type`" ``[string]`` When to use EWC on the predictor.  One
      of:
      
      - `"none`" Never do EWC
      - `"ewc`" Always do EWC
      - `"smart ewc`" Attempt EWC when it seems likely it will be useful and
        take the EWC correction if it is smaller than the standard correction.

    * `"freeze-thaw cusp width [K]`" ``[double]`` Controls a width over which
      to assume we are close to the latent heat cliff, and begins applying the
      EWC algorithm in `"ewc smarter`".

    * `"freeze-thaw cusp width (freezing) [K]`" ``[double]`` Controls a width
      over which to assume we are close to the latent heat cliff as we get
      colder, and begins applying the EWC algorithm in `"ewc smarter`".

    * `"freeze-thaw cusp width (thawing) [K]`" ``[double]`` Controls a width
      over which to assume we are close to the latent heat cliff as we get
      warmer, and begins applying the EWC algorithm in `"ewc smarter`".
        
    * `"pressure key`" ``[string]`` **DOMAIN-pressure**
    * `"temperature key`" ``[string]`` **DOMAIN-temperature**
    * `"water content key`" ``[string]`` **DOMAIN-water_content**
    * `"energy key`" ``[string]`` **DOMAIN-energy**
    * `"cell volume key`" ``[string]`` **DOMAIN-cell_volume**

    INCLUDES

    - ``[debugger-spec]`` Uses a Debugger_
    



State
#####
 State, a container for data.

State  is a  simple data-manager,  allowing PKs  to require,  read, and  write
various fields.

- Acts as a factory for data through the various require methods.
- Provides some data protection by providing both const and non-const
  data pointers to PKs.
- Provides some initialization capability -- this is where all
  independent variables can be initialized (as independent variables
  are owned by state, not by any PK).


.. _state-spec:
.. admonition:: state-spec

   * `"field evaluators`" ``[field-evaluator-typedinline-spec-list]`` A list of evaluators.
      Note this will eventually be an [evaluator-typedinline-spec-list] but the
      evaluators themselves do not include the type info.

   * `"initial conditions`" ``[list]`` A list of constant-in-time variables :
       `"initial conditions`" is a terrible name and will go away in the next
       iteration of state.

.. _field-evaluator-typedinline-spec:
.. admonition:: field-evaluator-typedinline-spec

   * `"field evaluator type`" ``[string]`` Type of the evaluator Included for
      convenience in defining data that is not in the dependency graph,
      constants are things (like gravity, or atmospheric pressure) which are
      stored in state but never change.  Typically they're limited to scalars
      and dense, local vectors.

.. _constants-scalar-spec:
.. admonition:: constants-scalar-spec

   * `"value`" ``[double]`` Value of a scalar constant

.. _constants-vector-spec:
.. admonition:: constants-vector-spec

   * `"value`" ``[Array(double)]`` Value of a dense, local vector.

Example:

.. code-block:: xml

    <ParameterList name="state">
      <ParameterList name="field evaluators">
        <ParameterList name="pressure">
          <Parameter name="field evaluator type" type="string" value="primary variable field evaluator" />
        </ParameterList>
      </ParameterList>

      <ParameterList name="initial conditions">
        <ParameterList name="atmospheric pressure">
          <Parameter name="value" type="double" value="101325.0" />
        </ParameterList>
        <ParameterList name="gravity">
          <Parameter name="value" type="Array(double)" value="{0.0,0.0,-9.80665}" />
        </ParameterList>
      </ParameterList>
    </ParameterList>




State consists of two sublists, one for evaluators and the other for
atomic constants.  The latter is currently called `"initial
conditions`", which is a terrible name which must be fixed.

example:

.. code-block:: xml
                
  <ParameterList name="state">
    <ParameterList name="field evaluators">
      ...
    </ParameterList>
    <ParameterList name="initial conditions">
      ...
    </ParameterList>
  </ParameterList>

 

Field Evaluators
=================
 A FieldEvaluator is a node in the dependency graph.

Evaluators combine data (variables) with functions (physics) providing the
values inside of that data.  Evaluators are the core concept of a model in
Amanzi-ATS: much if not all physics lives inside of an evaluator.

All evaluators are one of three types: primary, secondary, or independent
variable evaluators.  Primary variables are those that have no dependencies,
and are solved for.  Independent variables have no dependencies as well, but
are data provided by the user.  They may be constant or variable in time and
space, but they are assumed given.  Secondary variables are functions of time
and space, and of other variables.

From a software perspective, evaluators are nodes in a dependency graph, and
therefore provide the task ordering needed to evaluate the model.  Primary
variable evaluators do next to no work, but provide a target for solver to tell
the dependency graph when primary variables have been updated through the
solution/iteration of the nonlinear solve.  Independent and secondary variable
evaluators, on the other hand, implement the required processes needed to
update their data.

All field evaluators have keys they depend upon (unless they are a leaf) and
keys they provide.

All evaluator lists must provide an evaluator type, which is one of the types
registered with the evaluator factory.

.. _evaluator-typedinline-spec:
.. admonition:: evaluator-typedinline-spec

   * `"field evaluator type`" ``[string]`` Type registered in evaluator factory.
   * `"write checkpoint`" ``[bool]`` **true** Write this data when checkpointing.
   * `"write vis`" ``[bool]`` **true** Write this data when visualizing.





.. contents:: **List of Evalutors**
   :local:

PrimaryVariableEvaluator
------------------------


IndependentVariableEvaluator
----------------------------
 A field evaluator providing user-supplied data as a leaf node.

Independent variables are data provided by the user.  Independent variable
evaluators are therefore leaf nodes in the dependency graph (they have no
dependencies themselves) but do know how to update their data with values as a
function of space and time.

There are three ways for the user to provide that data: as static, constant
values, from an Amanzi Function_, or from file.




Constant
^^^^^^^^
 An evaluator with no dependencies specified by a constant value.

An independent variable evaluator that is constant in time.  This uses the same
infrastructure as initial conditions for time-independent data.

This is used by providing:

`"field evaluator type`" == `"independent variable constant`"

See the InitialConditions_ description for all parameters.



From Function
^^^^^^^^^^^^^
  A field evaluator with no dependencies specified by a function.

This evaluator is typically used for providing data that are functions of space
and time.  The evaluator consists of a list of region,function pairs, and the
functions are evaluated across that region at each timestep.  If the problem is
time-independent, the `"constant in time`" option results in a performance
boost (as the functions need only be evaluated once).  This leverages the
exaustive functional format capability provided in Amanzi's Functions_ library.

This evaluator is used by providing the option:

`"field evaluator type`" == `"independent variable`"

.. _independent-variable-evaluator-spec:
.. admonition:: independent-variable-evaluator-spec

   * `"constant in time`" ``[bool]`` **false** If true, only evaluate the
      functions once as they are time-independent.
   * `"function`" ``[composite-vector-function-spec-list]``

 


From File
^^^^^^^^^



SecondaryVariableEvaluators
---------------------------

All other evaluators are secondary variable evaluators, and these are grouped by physics concept or process type.

Conserved quantities
--------------------

Nearly all ATS process kernels are conservation equations, where there
is a formal conserved quantity, that, upon convergence, is conserved
to tolerance.  These are always an "extensive" quantity.

Water content in ATS is always measured on ``[mol]`` and therefore
includes a factor of the cell volume.  Energy in ATS is always
measured in ``[MJ]``.  Unlike nearly all other variables, this is not
SI, and is done so because this makes for fairly evenly balanced
equations between a coupled flow and energy problem.

Richards Equation water content (liquid only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Richards water content evaluator: the standard form as a function of liquid saturation.

.. math::
  \Theta = n s \phi V

Specified with evaluator type: `"richards water content`"

.. _field-evaluator-type-richards-water-content-spec:
.. admonition:: field-evaluator-type-richards-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"saturation liquid`"
   - `"cell volume`"




Liquid+Gas water content
^^^^^^^^^^^^^^^^^^^^^^^^
 Water content for liquid + water vapor.

.. math::
  \Theta = (n_l s_l + n_g s_g \omega) \phi V


Specified with evaluator type: `"liquid+gas water content`"

.. _field-evaluator-type-liquid-gas-water-content-spec:
.. admonition:: field-evaluator-type-liquid-gas-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"molar density gas`"
   - `"saturation liquid`"
   - `"saturation gas`"
   - `"mol frac gas`"
   - `"cell volume`"




Liquid+Ice water content
^^^^^^^^^^^^^^^^^^^^^^^^
 Water content for liquid + water vapor.

.. math::
  \Theta = (n_l s_l + n_i s_i) \phi V


Specified with evaluator type: `"liquid+ice water content`"

.. _field-evaluator-type-liquid-ice-water-content-spec:
.. admonition:: field-evaluator-type-liquid-ice-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"molar density ice`"
   - `"saturation liquid`"
   - `"saturation ice`"
   - `"cell volume`"




Liquid+Ice+Gas water content
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Three phase water content: vapor, liquid, and ice.

.. math::
  \Theta = (n_l s_l + n_i s_i + n_g s_g \omega_g ) \phi V

Specified with evaluator type: `"three phase water content`"

.. _field-evaluator-type-three-phase-water-content-spec:
.. admonition:: field-evaluator-type-three-phase-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"saturation liquid`"
   - `"molar density ice`"
   - `"saturation ice`"
   - `"molar density gas`"
   - `"saturation gas`"
   - `"molar fraction gas`"
   - `"cell volume`"




Surface water content
^^^^^^^^^^^^^^^^^^^^^


Snow or canopy water content
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most other water contents can be formed as `Multiplicative`_
evaluators.  See below for a few examples:

Multiplicative evalutator for `snow-water_content`:

.. code-block:: xml

      <ParameterList name="snow-water_content" type="ParameterList">
        <Parameter name="field evaluator type" type="string" value="multiplicative evaluator" />
        <Parameter name="evaluator dependencies" type="Array(string)" value="** DOC GENERATION ERROR: file not found 'snow-cell_volume, snow-water_equivalent, surface-molar_density_liquid' **" />
        <Parameter name="units" type="string" value="mol" />
      </ParameterList>

Multiplicative evaluator for `canopy-water_content`:

.. code-block:: xml

      <ParameterList name="canopy-water_content" type="ParameterList">
        <Parameter name="field evaluator type" type="string" value="multiplicative evaluator" />
        <Parameter name="evaluator dependencies" type="Array(string)" value="** DOC GENERATION ERROR: file not found 'canopy-cell_volume, canopy-water_equivalent, surface-molar_density_liquid' **" />
        <Parameter name="units" type="string" value="mol" />
      </ParameterList>

Richards energy
^^^^^^^^^^^^^^^
 Energy content evaluator for a standard Richards equation's water content.

Energy associated with a soil and pore-water system, in [KJ]

.. math::
  E = V * ( \phi u_l s_l n_l + (1-\phi_0) u_r \rho_r )

Specified with evaluator type: `"richards energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

Note that this ignores energy in the gas phase.

.. _field-evaluator-type-richards-energy-spec:
.. admonition:: field-evaluator-type-richards-energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"density rock`" Units may be either [kg m^-3] or [mol m^-3]
   - `"internal energy rock`" Units may be either [KJ kg^-1] or [KJ mol^-1],
     but must be consistent with the above density.
   - `"cell volume`" [m^3]




Liquid+Gas energy
^^^^^^^^^^^^^^^^^
 Energy content evaluator for a standard Richards equation, including energy in the gas phase.

Calculates energy, in [KJ], via the equation:

.. math::
  E = V * ( \phi (u_l s_l n_l + u_g s_g n_g)  + (1-\phi_0) u_r \rho_r )

Specified with evaluator type: `"liquid+gas energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

.. _field-evaluator-type-liquid-gas-energy-spec:
.. admonition:: field-evaluator-type-liquid-gas-energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density gas`" [mol m^-3]
   - `"saturation gas`" [-]
   - `"internal energy gas`" [KJ mol^-1]
   - `"density rock`" Units may be either [kg m^-3] or [mol m^-3]
   - `"internal energy rock`" Units may be either [KJ kg^-1] or [KJ mol^-1],
     but must be consistent with the above density.
   - `"cell volume`" [m^3]




Liquid+Ice energy
^^^^^^^^^^^^^^^^^
 Energy content evaluator for a two-phase system, including energy in an ice phase.

Calculates energy, in [KJ], via the equation:

.. math::
  E = V * ( \phi (u_l s_l n_l + u_i s_i n_i)  + (1-\phi_0) u_r \rho_r )

Specified with evaluator type: `"liquid+ice energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

Note that this ignores energy in the gas phase.

.. _field-evaluator-type-liquid-ice-energy-spec:
.. admonition:: field-evaluator-type-liquid-ice-energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density ice`" [mol m^-3]
   - `"saturation ice`" [-]
   - `"internal energy ice`" [KJ mol^-1]
   - `"density rock`" Units may be either [kg m^-3] or [mol m^-3]
   - `"internal energy rock`" Units may be either [KJ kg^-1] or [KJ mol^-1],
     but must be consistent with the above density.
   - `"cell volume`" [m^3]




Liquid+Ice+Gas energy
^^^^^^^^^^^^^^^^^^^^^
 Energy content evaluator for a three-phase, gas, liquid, ice system including the surrounding soil.

Calculates energy, in [KJ], via the equation:

.. math::
  E = V * ( \phi (u_l s_l n_l + u_i s_i n_i + u_g s_g n_g)  + (1-\phi_0) u_r \rho_r )

Specified with evaluator type: `"three phase energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

.. _field-evaluator-type-three-phase-energy-spec:
.. admonition:: field-evaluator-type-three-phase-energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density ice`" [mol m^-3]
   - `"saturation ice`" [-]
   - `"internal energy ice`" [KJ mol^-1]
   - `"molar density gas`" [mol m^-3]
   - `"saturation gas`" [-]
   - `"internal energy gas`" [KJ mol^-1]
   - `"density rock`" Units may be either [kg m^-3] or [mol m^-3]
   - `"internal energy rock`" Units may be either [KJ kg^-1] or [KJ mol^-1],
     but must be consistent with the above density.
   - `"cell volume`" [m^3]




Surface water+ice energy
^^^^^^^^^^^^^^^^^^^^^^^^
 Energy content for a surface water, partially frozen system.

The energy associated with ponded water, in [KJ], given by:

.. math::
  E = V * ( \eta h u_l n_l + (1 - \eta) h u_i n_i )

Specified with evaluator type: `"surface ice energy`"

.. _field-evaluator-type-surface-ice-energy-spec:
.. admonition:: field-evaluator-type-surface-ice-energy-spec

   DEPENDENCIES:

   - `"ponded depth`"  Height of water above the land surface [m]
   - `"unfrozen fraction`"  The fraction of unfrozen water ranges from 0 to 1. [-]
   - `"molar density liquid`" [mol m^-3]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density ice`" [mol m^-3]
   - `"internal energy ice`" [KJ mol^-1]
   - `"cell volume`" [m^2]





Subsurface flow evaluators
--------------------------

Assorted evaluators used for subsurface flow processes,
including water retention curves, compressible pore space, relative
permeability, and their frozen equivalents.

Many of these evaluators show up in nearly all ATS simulations, as
subsurface flow of water is the core process underlying all of ATS
physics.  For real examples, see `ats-demos <https://github.com/amanzi/ats-demos>`_

Capillary pressure
^^^^^^^^^^^^^^^^^^


Capillary pressure of liquid on ice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Water Retention Model and Relative Permeability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Evaluates saturation through water retention models.

Water Retention Models (WRMs) determine the saturation as a function of
pressure and the relative permeability as a function of saturation.  Most
commonly used in practice is the van Genuchten model, but others are available.

.. _wrm-evaluator-spec
.. admonition:: wrm-evaluator-spec

   * `"WRM parameters`" ``[WRM-typedinline-spec-list]``

   KEYS:
   - `"saturation`" **determined from evaluator name** The name
       of the liquid saturation -- typically this is determined from
       the evaluator name and need not be set.
   - `"other saturation`"  **determined from evaluator name**
       The name of the other saturation, usually gas -- typically this is determined
       from the evaluator name and need not be set.
   - `"capillary pressure`"` **DOMAIN-capillary_pressure_gas_liq**
       The name of the capillary pressure.



 A collection of WRMs along with a Mesh Partition.

A WRM partition is a list of (region, WRM) pairs, where the regions partition
the mesh.

.. _wrm-partition-typedinline-spec
.. admonition:: wrm-partition-typedinline-spec

   * `"region`" ``[string]`` Region on which the WRM is valid.
   * `"WRM type`" ``[string]`` Name of the WRM type.
   * `"_WRM_type_ parameters`" ``[_WRM_type_-spec]`` See below for the required
      parameter spec for each type.



 Evaluates relative permeability using water retention models.

Uses a list of regions and water retention models on those regions to evaluate
relative permeability, typically as a function of liquid saturation.

Most of the parameters are provided to the WRM model, and not the evaluator.
Typically these share lists to ensure the same water retention curves, and this
one is updated with the parameters of the WRM evaluator.  This is handled by
flow PKs.

Some additional parameters are available.

.. _rel-perm-evaluator-spec
.. admonition:: rel-perm-evaluator-spec

   * `"use density on viscosity in rel perm`" ``[bool]`` **true**

   * `"boundary rel perm strategy`" ``[string]`` **boundary pressure** Controls
     how the rel perm is calculated on boundary faces.  Note, this may be
     overwritten by upwinding later!  One of:

      - `"boundary pressure`" Evaluates kr of pressure on the boundary face, upwinds normally.
      - `"interior pressure`" Evaluates kr of the pressure on the interior cell (bad idea).
      - `"harmonic mean`" Takes the harmonic mean of kr on the boundary face and kr on the interior cell.
      - `"arithmetic mean`" Takes the arithmetic mean of kr on the boundary face and kr on the interior cell.
      - `"one`" Sets the boundary kr to 1.
      - `"surface rel perm`" Looks for a field on the surface mesh and uses that.

   * `"minimum rel perm cutoff`" ``[double]`` **0.** Provides a lower bound on rel perm.

   * `"permeability rescaling`" ``[double]`` Typically rho * kr / mu is very big
     and K_sat is very small.  To avoid roundoff propagation issues, rescaling
     this quantity by offsetting and equal values is encourage.  Typically 10^7 or so is good.

   * `"WRM parameters`" ``[wrm-typedinline-spec-list]``  List (by region) of WRM specs.

   KEYS:

   - `"rel perm`"
   - `"saturation`"
   - `"density`" (if `"use density on viscosity in rel perm`" == true)
   - `"viscosity`" (if `"use density on viscosity in rel perm`" == true)
   - `"surface relative permeability`" (if `"boundary rel perm strategy`" == `"surface rel perm`")




Van Genuchten Model
~~~~~~~~~~~~~~~~~~~
 WRMVanGenuchten : water retention model using van Genuchten's parameterization

van Genuchten's water retention curve.

.. _WRM-van-Genuchten-spec
.. admonition:: WRM-van-Genuchten-spec

    * `"region`" ``[string]`` Region to which this applies
    * `"van Genuchten alpha [Pa^-1]`" ``[double]`` van Genuchten's alpha

    ONE OF:
    * `"van Genuchten n [-]`" ``[double]`` van Genuchten's n
    OR
    * `"van Genuchten m [-]`" ``[double]`` van Genuchten's m, m = 1 - 1/n
    END

    * `"residual saturation [-]`" ``[double]`` **0.0**
    * `"smoothing interval width [saturation]`" ``[double]`` **0.0**
    * `"Mualem exponent l [-]`" ``[double]`` **0.5**
    * `"Krel function name`" ``[string]`` **Mualem**  `"Mualem`" or `"Burdine`"

Example:

.. code-block:: xml

    <ParameterList name="moss" type="ParameterList">
      <Parameter name="region" type="string" value="moss" />
      <Parameter name="WRM type" type="string" value="van Genuchten" />
      <Parameter name="van Genuchten alpha [Pa^-1]" type="double" value="0.002" />
      <Parameter name="van Genuchten m [-]" type="double" value="0.2" />
      <Parameter name="residual saturation [-]" type="double" value="0.0" />
      <Parameter name="smoothing interval width [saturation]" type="double" value=".05" />
    </ParameterList>




Linear  Model
~~~~~~~~~~~~~~~~~~~
 A linear sat-pc curve.

  A linear sat-pc curve, plus a constant rel perm, makes the system linear, so
  nonlinear solver should always converge in one step.

  No error-checking, so the user is responsible for ensuring that the pressure
  is always less than atmospheric and within the acceptable range of the slope.

  Note this is mostly for testing.




Water Retention Model for Freeze-Thaw
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Original Implicit model
~~~~~~~~~~~~~~~~~~~~~~~
 Painter's original, implicitly defined permafrost model.

.. _wrm-implicit-permafrost-spec
.. admonition:: wrm-implicit-permafrost-spec

    * `"converged tolerance`" ``[double]`` **1.e-12** Convergence tolerance of the implicit solve.
    * `"max iterations`" ``[int]`` **100** Maximum allowable iterations of the implicit solve.
    * `"solver algorithm [bisection/toms]`" ``[string]`` **bisection** Use bisection or the TOMS algorithm from boost.




Freezing point depression model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Freezing point depression, smoothed model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** DOC GENERATION ERROR: file not found ' wrm_fpd_smoothed_model ' **

Interfrost model
~~~~~~~~~~~~~~~~


Sutra-ICE model
~~~~~~~~~~~~~~~


Compressible porosity
^^^^^^^^^^^^^^^^^^^^^


Compressible grains are both physically realistic (based on bulk modulus) and a
simple way to provide a non-elliptic, diagonal term for helping solvers to
converge.

* `"compressible porosity model parameters`" ``[compressible-porosity-model-spec-list]``

KEYS:
- `"pressure`" **DOMAIN-pressure**
- `"base porosity`" **DOMAIN-base_porosity**




Standard model
~~~~~~~~~~~~~~
 A simple model for allowing porosity to vary with pressure.

Based on a linear increase, i.e.

$$ \phi = \phi_{base} + H(p - p_{atm}) * \alpha $$

where $H$ is the heaviside function and $\alpha$ is the provided
compressibility.  If the inflection point is set to zero, the above function
is exact.  However, then the porosity function is not smooth (has
discontinuous derivatives).

* `"pore compressibility [Pa^-1]`" ``[double]``  $\alpha$ as described above
* `"pore compressibility inflection point [Pa]`" ``[double]`` **1000**
* `"region`" ``[string]`` Region on which this is applied.

  The inflection point above which the function is linear.

Example:

.. code-block::xml

  <ParameterList name="soil" type="ParameterList">
    <Parameter name="region" type="string" value="soil" />
    <Parameter name="pore compressibility [Pa^-1]" type="double" value="1.e-9" />
    <Parameter name="pore compressibility inflection point [Pa]" type="double" value="1000." />
  </ParameterList>




Exponential model
~~~~~~~~~~~~~~~~~


The Leinjnse model is an exponential model of porosity as a function of
pressure, based on (insert citation!):

.. math:
   p_\lambda = p - p_\text{ref}
   \phi = 1 - (1-\phi_\text{base}) * (exp(- \alpha (p_\lambda - \delta)) - 0.5 \alpha \delta), p_\lambda > \delta
   \phi = 1 - (1-\phi_\text{base}) * (1 - 0.5 \alpha / \delta p_\lambda^2)

where :math:`\alpha` is the provided
compressibility, and :math:`\delta` is the cutoff (inflection point). 

If the inflection point is set to zero, the above function is exact.  However,
then the porosity function is not smooth (has discontinuous derivatives).
  
* `"pore compressibility [Pa^-1]`" ``[double]`` :math:`\alpha` as described above
  
* `"pore compressibility inflection point [Pa]`" ``[double]`` **1000** The inflection point above which the function is linear.

NOTE: provide a parameter in the `EWC Globalization Delegate`_ to turn Leijnse model ON 

  <Parameter name="porosity leijnse model" type="bool" value="true"/>





Viscosity Constant
^^^^^^^^^^^^^^^^^^
Like any quantity, a viscosity can simply be a constant value, at which
point it is not a secondary variable but an independent variable.

.. code-block:: xml

      <ParameterList name="viscosity_liquid" type="ParameterList">
        <Parameter name="field evaluator type" type="string" value="independent variable constant" />
        <Parameter name="value" type="double" value="8.9e-4" />
        <Parameter name="units" type="string" value="Pa s" />
      </ParameterList>

Viscosity of Water
^^^^^^^^^^^^^^^^^^




Surface flow evaluators
-----------------------

Assorted evaluators used for surface flow, including potential
surfaces, Manning's conductivity, and their frozen equivalents.

Like the subsurface flow evaluators, many of these evaluators show up
in nearly all ATS simulations.  For real examples, see `ats-demos
<https://github.com/amanzi/ats-demos>`_


Ponded Depth or Water Height
^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Ponded Depth, Frozen
^^^^^^^^^^^^^^^^^^^^



Effective, or Smoothed Height
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Unfrozen fraction
^^^^^^^^^^^^^^^^^



Unfrozen Flowing Depth
^^^^^^^^^^^^^^^^^^^^^^



SurfacePotential
^^^^^^^^^^^^^^^^^^^
 PresElevEvaluator: evaluates h + z
Evaluator type: ""

.. math::
  h + z

* `"my key`" ``[string]`` **pres_elev** Names the surface water potential variable, h + z [m]
* `"height key`" ``[string]`` **ponded_depth** Names the height variable. [m]
* `"elevation key`" ``[string]`` **elevation** Names the elevation variable. [m]


NOTE: This is a legacy evaluator, and is not in the factory, so need not be in
the input spec.  However, we include it here because this could easily be
abstracted for new potential surfaces, kinematic wave, etc, at which point it
would need to be added to the factory and the input spec.

NOTE: This could easily be replaced by a generic Additive_ Evaluator.




Overland Conductivity, sheet flow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. _`Overland Conductivity Evaluator`:



Overland Conductivity, litter resistance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




Thermodynamic evaluators
-------------------------

Internal energy, enthalpy, thermal conductivity, etc used for both
surface and subsurface transport of energy.

Internal energy
^^^^^^^^^^^^^^^


Linear
~~~~~~
 Internal energy based on a linear fit.

Linear internal energy model -- function of Cv and temperature

.. math::

    u = L_f +  C_v * (T - T_{ref})

.. _iem-linear-spec
.. admonition:: iem-linear-spec

    * `"reference temperature [K]`" ``[double]`` **273.15** The phase transition point, T_ref above

    ONE OF

    * `"latent heat [J kg^-1]`" ``[double]`` Latent heat of fusion, L_f above
    * `"heat capacity [J kg^-1 K^-1]`" ``[double]`` C_v above

    OR

    * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of fusion, L_f above.
    * `"heat capacity [J mol^-1 K^-1]`" ``[double]`` C_v above

    END




Quadratic
~~~~~~~~~
 Internal energy based on a quadratic fit to data.

Quadratic internal energy model -- function of Cv and temperature

.. math::

    u = L_f + C_v * (T - T_{ref}) + b(T - T_{ref})^2

.. _iem-quadratic-spec
.. admonition:: iem-quadratic-spec

    * `"reference temperature [K]`" ``[double]`` **273.15** The phase transition point, T_ref above.

    ONE OF

    * `"latent heat [J kg^-1]`" ``[double]`` Latent heat of fusion, L_f above
    * `"heat capacity [J kg^-1 K^-1]`" ``[double]`` C_v above
    * `"quadratic b [J kg^-1 K^-2]`" ``[double]`` b above

    OR

    * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of fusion, L_f above.
    * `"heat capacity [J mol^-1 K^-1]`" ``[double]`` C_v above
    * `"quadratic b [J mol^-1 K^-2]`" ``[double]`` b above

    END




Water Vapor
~~~~~~~~~~~
 Internal energy model for air and water vapor.

Internal energy model for air and water vapor, relative to water @237.15K

.. _iem-water-vapor-spec
.. admonition:: iem-water-vapor-spec

    * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of vaporization
    * `"heat capacity [J mol^-1 K^-1]`" ``[double]`` C_v




Enthalpy
~~~~~~~~


Thermal Conductivity, two phases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Wet-Dry Model
~~~~~~~~~~~~~


Peters-Lidard Model
~~~~~~~~~~~~~~~~~~~


Thermal Conductivity, three phases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Thermal conductivity based on a three-phase (air,liquid,ice) composition of the porous media.

Wet-Dry Model
~~~~~~~~~~~~~
 Three-phase thermal conductivity based on paper by Peters-Lidard.

A three-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Empirical relationship for frozen soil based on Peters-Lidard

See Atchley et al GMD 2015 Supplementary Material for equations.

.. _thermal-conductivity-threephase-wetdry-spec:
.. admonition:: thermal-conductivity-threephase-wetdry-spec

    * `"thermal conductivity, saturated (unfrozen) [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of fully saturated, unfrozen bulk soil.
    * `"thermal conductivity, dry [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of fully dried bulk soil.
    * `"unsaturated alpha unfrozen [-]`" ``[double]`` Interpolating exponent
    * `"unsaturated alpha frozen [-]`" ``[double]`` Interpolating exponent
    * `"unsaturated alpha frozen [-]`" ``[double]`` Interpolating exponent
    * `"saturated beta frozen [-]`" ``[double]`` **1.0** Interpolating exponent
    * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.




Peters-Lidard Model
~~~~~~~~~~~~~~~~~~~
 Three-phase thermal conductivity based on paper by Peters-Lidard.

A three-phase thermal conductivity, based upon:

- A mixture model using interpolation across various components.
- Power-law Kersten number.

See Atchley et al GMD 2015 Supplementary Material for equations.

.. _thermal-conductivity-threephase-peterslidard-spec:
.. admonition:: thermal-conductivity-threephase-peterslidard-spec

    * `"thermal conductivity of soil [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of soil grains (not bulk soil)
    * `"thermal conductivity of liquid [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of liquid (water)
    * `"thermal conductivity of gas [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of gas (air)
    * `"thermal conductivity of ice [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of ice
    * `"unsaturated alpha unfrozen [-]`" ``[double]`` Interpolating exponent
    * `"unsaturated alpha frozen [-]`" ``[double]`` Interpolating exponent
    * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

Usage:

  <ParameterList name="thermal_conductivity">
    <Parameter name="thermal conductivity type" type="string" value="three-phase Peters-Lidard"/>
    <Parameter name="thermal conductivity of soil [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of gas [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of ice [W m^-1 K^-1]" type="double" value=""/>

    <Parameter name="unsaturated alpha unfrozen [-]" type="double" value=""/>
    <Parameter name="unsaturated alpha frozen [-]" type="double" value=""/>

    <Parameter name="epsilon" type="double" value="1.e-10"/>
  </ParameterList>




Volume-averaged Model
~~~~~~~~~~~~~~~~~~~~~
 A volume-averaged thermal conductivity based on TCs of raw components.

A simple model of three-phase thermal conductivity, based upon volume-averaging
of four consitutive components.

See Atchley et al GMD 2015 Supplementary Material for equations.

.. _thermal-conductivity-volume-averaged-spec:
.. admonition:: thermal-conductivity-volume-averaged-spec

    * `"thermal conductivity of soil [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of soil **grains**
    * `"thermal conductivity of liquid [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of liquid water.
    * `"thermal conductivity of gas [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of air.
    * `"thermal conductivity of ice [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of frozen water.




Sutra-ICE model
~~~~~~~~~~~~~~~~~~~


Thermal Conductivity, Surface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Advected Energy Source
^^^^^^^^^^^^^^^^^^^^^^


Active Layer Averaged Temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Water Table
^^^^^^^^^^^
 An evaluator for calculating the depth to water table, relative to the surface.

Note that postive values are below the surface, and negative values
are not possible.  One may choose to, in post-processing, do:

  water_table_height = where(ponded_depth > 0, ponded_depth, -water_table_depth)

to calculate a water table height that is valid for both surface and
subsurface water.

Evaluator name: `"water table depth`"

.. _water-table-depth-spec:
.. admonition:: water-table-depth-spec

    KEYS:
      `"saturation_gas`"




Thaw Depth
^^^^^^^^^^
 An evaluator for calculating the depth to frozen soil/permafrost, relative to the surface.

Note that postive values are below the surface, and negative values
are not possible. 

Evaluator name: `"thaw depth`"

.. _thaw-depth-spec:
.. admonition:: thaw-depth-spec

    * `"transition width [K]`" ``[double]`` Width of the freeze curtain transition.

    KEYS:
      `"temperature`"





Equations of State
------------------

The density of water can be specified in many ways, depending upon
phase and problem of interest.  Options are available from the
simplest (constant value) to functions of just temperature,
temperature and pressure, and temperature/pressure/concentration
(e.g. salinity).

Note that density includes both molar and mass-based values.  Most
density evaluators can provide either, or do provide both, the
difference being simply a factor of the molecular mass of water.

Finally, nearly all (except the constant value) equations of state use
a common set of models which accept all of temperature, pressure, and
concentration.  Many of these models ignore one or the other, so it
should be preferred (but is not necessary) to choose the right model
for a given evaluator.  Choosing a model that uses fewer of
these than the evaluator provides is valid, but it is inefficient.
Choosing a model that uses more than the evaluator provides will
result in an error.

Constant Value
^^^^^^^^^^^^^^

Like any quantity, a density can simply be a constant value, at which
point it is not a secondary variable but an independent variable.

.. code-block:: xml

      <ParameterList name="surface-molar_density_liquid" type="ParameterList">
        <Parameter name="field evaluator type" type="string" value="independent variable constant" />
        <Parameter name="value" type="double" value="55000" />
        <Parameter name="units" type="string" value="mol m^-3" />
      </ParameterList>

Isobaric (constant pressure)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Standard Equation of State
^^^^^^^^^^^^^^^^^^^^^^^^^^


Conc-Temp-Pres Equation of State
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Models
^^^^^^

Constant EOS
~~~~~~~~~~~~


Linear EOS
~~~~~~~~~~


Ideal Gas
~~~~~~~~~


EOS of Water
~~~~~~~~~~~~


EOS of Ice
~~~~~~~~~~


EOS of Vapor in Air
~~~~~~~~~~~~~~~~~~~


EOS of Saltwater
~~~~~~~~~~~~~~~~


Surface energy balance evaluators
---------------------------------

Evaluators used to solve the fluxes to and from the atmosphere and
between layers of the surface.  Typically in ATS these calculate
evapotranspiration.

Area Fractions
^^^^^^^^^^^^^^

Frequently, the surface of a grid cell is split across at least two
"regions," for instance snow covered and bare ground.  This "subgrid"
model allows smooth transitions between fully snow-covered and fully
bare ground.

These area fractions often get included as area-weights in calculating
full-cell quantities.

Surface Area Fractions
~~~~~~~~~~~~~~~~~~~~~~
** DOC GENERATION ERROR: file not found ' area_fractions_evaluator ' **


Potential Evapotranspiration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Models of potential evapotranspiration approximate the difference in
vapor pressure between the atmosphere and the soil as a function of
available energy, allowing the calculation of the max flux of ET that
the atmosphere could accept.  This can then be limited based on water
availability, etc.

Priestley-Taylor Potential Evapotranspiration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** DOC GENERATION ERROR: file not found ' pet_evaluator ' **

Downregulation and limiters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a potential, the actual ET is often limited by available water
(or nutrients or other quantities).  These evaluators are used to
limit, downregulate, distribute, or otherwise move a potential to an
actual ET.

Transpiration Distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~
 Distributes and downregulates potential transpiration to the rooting zone.

The transpiration distribution evaluator looks to take a potential
evapotranspiration and distribute it across the vertical column based on water
availability and rooting depths.  It also potentially limits the transpiration
to avoid taking water where it is not available (thereby crashing the code).

This model is based off of versions from both CLM 4.5 and PRMS.  It requires:

1. A root distribution profile.
2. A plant wilting factor (e.g. how water stressed is the plant?)
3. A potential transpiration, typically calculated from data or a potential
   difference based on a latent heat calculation.

Note this also requires columnar meshes -- meaning that the subsurface mesh
must have `"build columns from set`" provided.

A normalized fraction of water is calculated through multiplying the water
factor by the root distribution factor, integrating across the column, and
dividing by the integral.  This gives a factor which sums to 1 and can be used
to distribute the potential ET throughout the soil column.

Then, this potential ET is down-regulated and multiplied by the plant wilting
factor.  If there is no water locally, it cannot be taken.  Note that almost
always, if there is no water, this did not contribute (much) to the integral
and so is already small.  But if the entire root zone is dried out, this might
have been a bunch of small numbers which then got normalized to 1, meaning they
are all now significant.

Finally, transpiration may be turned off for winter -- relative to time zero,
parameters `"leaf on doy`" and `"leaf off doy`" are used to control when ET
is zero.  By default these are set to 0 and 366 days, ensuring that
transpiration is never turned off and the potential ET must include this
factor.  This is the right choice for, e.g. ELM output, or eddy covariance flux
tower data (where leaf on and off are already included in the potential
calculation).  It is the wrong choice for, e.g. Priestly-Taylor or
Penmann-Montief models, which are happy to predict transpiration all winter
long.  Good choices for those models depend upon the local climate, but may be
something like Julian day 101 for leaf on and Julian day 254 for leaf off (PRMS
defaults for US temperate forests).

Note that `"leaf on doy`" and `"leaf off doy`" are relative to the simulation's
zero time, not the start time.  Typically these are Julian day of the year, but
this assumes that the 0 time of the simulation (not the "start time", but time
0!) is Jan 1.  This leaf on/off cycle is modulo the `"year duration`"
(typically 1 noleap).  Note if `"leaf off doy`" < `"leaf on time`" is ok too --
this is the case if simulation time zero is mid-summer.  These parameters come
from the LandCover type.

.. _transpiration-distribution-evaluator-spec:
.. admonition:: transpiration-distribution-evaluator-spec

    * `"year duration`" ``[double]`` **1**
    * `"year duration units`" ``[string]`` **noleap**

    * `"water limiter function`" ``[function-spec]`` **optional** If provided,
      limit the total water sink as a function of the integral of the water
      potential * rooting fraction.

    KEYS:

    - `"plant wilting factor`" **DOMAIN-plant_wilting_factor**
    - `"rooting depth fraction`" **DOMAIN-rooting_depth_fraction**
    - `"potential transpiration`" **DOMAIN_SURF-potential_transpiration**
    - `"cell volume`" **DOMAIN-cell_volume**
    - `"surface cell volume`" **DOMAIN_SURF-cell_volume**




Rooting Depth Fraction
~~~~~~~~~~~~~~~~~~~~~~
 Provides a depth-based profile of root density.

Sets the root fraction as a function of depth,

.. math:
   F_root =  ( \alpha \; exp(-\alpha z) + \beta \; exp(-\beta z) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

Note that all three parameters, a, b, and the cutoff, are provided in the
LandCover type.

.. _rooting-depth-fraction-evaluator-spec:
.. admonition:: rooting-depth-fraction-evaluator-spec

   * `"surface domain name`" ``[string]`` **SURFACE_DOMAIN** Sane default provided for most domain names.

   KEYS:
   - `"depth`" **DOMAIN-depth**
   - `"cell volume`" **DOMAIN-cell_volume**
   - `"surface cell volume`" **SURFACE_DOMAIN-cell_volume**




Plant Wilting Point
~~~~~~~~~~~~~~~~~~~
 Plant wilting factor provides a moisture availability-based limiter on transpiration.

  Generated via evaluator_generator with:
Wilting factor.

Beta, or the water availability factor, or the plant wilting factor.

Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or soil mafic potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).

Note this makes use of LandCover objects for mafic potential of fully open and
fully closed stomata.

Note the challenges of using this model with arbitrary van Genuchten WRMs.  See
Verhoef & Egea, Ag. & Forest Meteorology, 2014
https://doi.org/10.1016/j.agrformet.2014.02.009


.. _plant-wilting-factor-evaluator-spec:
.. admonition:: plant-wilting-factor-evaluator-spec

   KEYS:
   - `"capillary pressure`" **DOMAIN-capillary_pressure_gas_liq**





Soil Resistance
^^^^^^^^^^^^^^^
 Downregulates evaporation via vapor diffusion through a dessicated zone.

Calculates evaporative resistance through a dessicated zone.

Sakagucki and Zeng 2009 equations 9 and 10.

Requires the use of LandCover types, for dessicated zone thickness and Clapp &
Hornberger b.

.. _evaporation-downregulation-evaluator-spec:
.. admonition:: evaporation-downregulation-evaluator-spec

   KEYS:
   - `"saturation gas`" **DOMAIN_SUB-saturation_gas**
   - `"porosity`" **DOMAIN_SUB-porosity**
   - `"potential evaporation`" **DOMAIN_SUB-potential_evaporation**





Radiation Balance Terms
^^^^^^^^^^^^^^^^^^^^^^^

Often a balance of incoming and outgoing short and longwave radiations
are required to determine the energy available to go into latent heat,
and therefore potential evapotranspiration.

Note that incoming shortwave radiation is typically a user-provided
meterological forcing dataset.

Surface Albedo
~~~~~~~~~~~~~~
** DOC GENERATION ERROR: file not found ' albedo_evaluator ' **

Incident Shortwave Radiation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Evaluates shortwave as a function of slope/aspect/etc.
  Generated via evaluator_generator with:
Aspect modified shortwave radiation is determined by a factor which
is multiplied by the 'incoming radiation incident on a flat surface'
to determine the 'incoming radiation incident on a sloping surface of
a given aspect' as a function of latitude, slope, aspect, and Julian
day of the year, and time of day.

Note that some careful checking and experimentation has found that, in
general, the daily average incident radiation times the 12-noon aspect
modifier correlates reasonably well with the daily average of the
product of the hourly incident radiation and the hourly aspect
modifier.  It is notably better than the daily average radiation times
the daily average aspect modifier.

Derived from LandLab code, which is released under the MIT license:
https://github.com/landlab/landlab/blob/master/landlab/components/radiation/radiation.py



Longwave Radiation
~~~~~~~~~~~~~~~~~~
 Evaluates incoming longwave radiation from rel humidity and air temperature.

.. _longwave_evaluator-spec:
.. admonition:: longwave_evaluator-spec

    * `"minimum relative humidity [-]`" ``[double]`` **0.1** Sets a minimum rel humidity, RH=0 breaks the model.

    DEPENDENCIES:

    * `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
    * `"relative humidity key`" ``[string]`` **DOMAIN-relative_humidity**




Full Surface Energy Balance Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, in addition to the potential-based models above, a few
full-physics model are available.  These are often based on older,
monolithic models.

Bare Soil Surface Energy Balance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** DOC GENERATION ERROR: file not found ' seb_evaluator ' **

Common Land Model (ParFlow-CLM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is included here because it should eventually get split into
evaluators.  Currently, the CLM SEB model is a PK only, see `Common
Land Model PK`_.

Snow evaluators
---------------

Evaluators used for solving snowpack and/or snow redistribution laterally.

SnowSurfacePotential
^^^^^^^^^^^^^^^^^^^^^^
 PresElevEvaluator: evaluates h + z
Evaluator type: "snow skin potential"

.. math::
  h + z + h_{{snow}} + dt * P_{{snow}}

* `"my key`" ``[string]`` **snow_skin_potential** Names the potential variable evaluated [m]
* `"ponded depth key`" ``[string]`` **ponded_depth** Names the surface water depth variable. [m]
* `"snow depth key`" ``[string]`` **snow_depth** Names the snow depth variable. [m]
* `"precipitation snow key`" ``[string]`` **precipitation_snow** Names the snow precipitation key. [m]
* `"elevation key`" ``[string]`` **elevation** Names the elevation variable. [m]
* `"dt factor [s]`" ``[double]`` A free-parameter factor for providing a time scale for diffusion of snow precipitation into low-lying areas.  Typically on the order of 1e4-1e7. This timestep times the wave speed of snow provides an approximate length of how far snow precip can travel.  Extremely tunable! [s]

NOTE: This is equivalent to a generic Additive_ Evaluator

Example:

.. code-block:: xml

  <ParameterList name="snow_skin_potential" type="ParameterList">
    <Parameter name="field evaluator type" type="string" value="snow skin potential" />
    <Parameter name="dt factor [s]" type="double" value="864000.0" />
  </ParameterList>




Snow Melt Rate
^^^^^^^^^^^^^^
 Evaluates snow melt via USDA - Natural Resources Conservation Service model

From:  National Engineering Handbook (NEH) part 630, Chapter 11

Uses LandCover for snow_ground_transition parameter.

.. _snow-meltrate-evaluator-spec:
.. admonition:: snow-meltrate-evaluator-spec

   * `"snow melt rate [mm day^-1 C^-1]`" ``[double]`` **2.74**
     the melt rate per degree-day above 0 C.

   * `"air-snow temperature difference [C]`" ``[double]`` **2.0**
     Snow temperature is typicaly a few degrees colder than the air
     temperature at snowmelt. This shifts air temp (positive is colder)
     when calculating the snow temperature.

   * `"surface domain name`" ``[string]`` **SURFACE_DOMAIN** Attempts to
     guess a sane default by the snow domain name.

   KEYS:
   - `"air temperature`"  **SURFACE_DOMAIN-air_temperature**
   - `"snow water equivalent`" **DOMAIN-water_equivalent**


.. note:
    If snow temperature is known, the `"air-snow temperature difference`"
    should be set to 0, and the `"air temperature key`" should be set to
    the snow temperature key instead.





Subgrid flow model evaluators
-----------------------------

Surface flow is sensitive to microtopographic effects, and the surface
flow model of Jan2018_ can be used to represent these effects.

Depression Depth
^^^^^^^^^^^^^^^^


Excluded Volume
^^^^^^^^^^^^^^^


Microtopographic Relief
^^^^^^^^^^^^^^^^^^^^^^^


Drag Exponent
^^^^^^^^^^^^^


Surface Albedo with subgrid
^^^^^^^^^^^^^^^^^^^^^^^^^^^
** DOC GENERATION ERROR: file not found ' albedo_subgrid_evaluator ' **

Surface Area Fractions with subgrid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
** DOC GENERATION ERROR: file not found ' area_fractions_subgrid_evaluator ' **


Deformation
-----------
Evaluators for implementing ATS's simple deformation model targetting subsidence.

Initial Elevation
^^^^^^^^^^^^^^^^^


Subsidence
^^^^^^^^^^


Deforming Porosity
^^^^^^^^^^^^^^^^^^



Carbon and Biogeochemistry Models
---------------------------------

Carbon Decomposition Rate
^^^^^^^^^^^^^^^^^^^^^^^^^


Moisture Content
^^^^^^^^^^^^^^^^


Cryo/bioturbation
^^^^^^^^^^^^^^^^^


Pool Decomposition
^^^^^^^^^^^^^^^^^^


Pool Transfer
^^^^^^^^^^^^^



Multiscale Models
-----------------

Subgrid Aggregation
^^^^^^^^^^^^^^^^^^^^^^
 SubgridAggregateEvaluator restricts a field to the subgrid version of the same field.

.. _subgrid-aggregate-evaluator-spec:
.. admonition:: subgrid-aggregate-evaluator-spec

   * `"source domain name`" ``[string]`` Domain name of the source mesh.

   KEYS:
   - `"field`" **SOURCE_DOMAIN-KEY**  Default set from this evaluator's name.




Subgrid disaggregation
^^^^^^^^^^^^^^^^^^^^^^
 SubgridDisaggregateEvaluator restricts a field to the subgrid version of the same field.

.. _subgrid-disaggregate-evaluator-spec:
.. admonition:: subgrid-disaggregate-evaluator-spec

   * `"source domain name`" ``[string]`` Domain name of the source mesh.

   KEYS:
   - `"field`" **SOURCE_DOMAIN-KEY**  Default set from this evaluator's name.


 






Geometric evaluators
--------------------

Evaluators that capture various aspects of the mesh geometry.

SurfaceElevation
^^^^^^^^^^^^^^^^^^
 MeshedElevationEvaluator: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.
Evaluator type: `"meshed elevation`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

* `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
* `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation variable. [-]
* `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the elevation changes in time, and adds the `"deformation`" dependency.
* `"parent domain name`" ``[string]`` **DOMAIN** Domain name of the parent mesh, which is the 3D version of this domain.  Attempts to generate an intelligent default by stripping "surface" from this domain.

Example:

.. code-block:: xml

  <ParameterList name="elevation">
    <Parameter name="evaluator type" type="string" value="meshed elevation"/>
  </ParameterList>




Elevation of a Column
^^^^^^^^^^^^^^^^^^^^^
 ElevationEvaluatorColumn: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.
Evaluator type: `"elevation column`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

* `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
* `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation variable. [-]
* `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the elevation changes in time, and adds the `"deformation`" and `"base_porosity`" dependencies.
* `"parent domain name`" ``[string]`` **DOMAIN** Domain name of the parent mesh, which is the 3D version of this domain.  In the columnar meshes the surface elevation and slope are assigned based on the columns and not the base 3D domain.

Example:

.. code-block:: xml

  <ParameterList name="surface_star-elevation">
    <Parameter name="field evaluator type" type="string" value="column elevation"/>
  </ParameterList>




Depth Evaluator
^^^^^^^^^^^^^^^


Cell Volume evaluator
^^^^^^^^^^^^^^^^^^^^^



Coupling terms
--------------


Generic evaluators
------------------

These evaluators take arguments to allow the user to provide custom
functions via the input spec.

Additive
^^^^^^^^
 A generic evaluator for summing a collection of fields.

.. _additive-evaluator-spec:
.. admonition:: additive-evaluator-spec
   * `"constant shift`" ``[double]`` **0** A constant value to add to the sum.

   * `"DEPENDENCY coefficient`" ``[double]`` A multiple for each dependency in
     the list below.

   ONE OF
   * `"evaluator dependencies`" ``[Array(string)]`` The things to sum.
   OR
   * `"evaluator dependency suffixes`" ``[Array(string)]``
   END




Multiplicative
^^^^^^^^^^^^^^
 A generic evaluator for multiplying a collection of fields.

.. _multiplicative-evaluator-spec:
.. admonition:: multiplicative-evaluator-spec
   * `"coefficient`" ``[double]`` **1** A constant prefix to the product.
   * `"enforce positivity`" ``[bool]`` **false** If true, max the result with 0.

   ONE OF
   * `"evaluator dependencies`" ``[Array(string)]`` The fields to multiply.
   OR
   * `"evaluator dependency suffixes`" ``[Array(string)]``
   END




Column summation
^^^^^^^^^^^^^^^^
 Sums a subsurface field vertically only a surface field.

Simple vertical sum of all cells below each surface cell.  Note that their are
options for including volume factors (multiply by subsurface cell volume, sum,
divide by surface cell area) and density (useful for the most common use case
of summing fluxes onto the surface and converting to m/s instead of mol/m^2/s).


.. _column-sum-evaluator-spec:
.. admonition:: column-sum-evaluator-spec

    * `"include volume factor`" ``[bool]`` **true** In summing, multiply the
      summand subsurface cell volume, then divide the sum by the surface cell
      area.  Useful for converting volumetric fluxes to total fluxes.

    * `"divide by density`" ``[bool]`` **true** Divide the summand by density.
      Useful for converting molar fluxes to volumetric fluxes
      (e.g. transpiration).

    * `"column domain name`" ``[string]`` **domain** The domain of the
      subsurface mesh.  Note this defaults to a sane thing based on the
      variable's domain (typically "surface" or "surface_column:*") and is
      rarely set by the user.

    KEYS:
    - `"summed`" The summand, defaults to the root suffix of the calculated variable.
    - `"cell volume`" Defaults to domain's cell volume.
    - `"surface cell volume`" Defaults to surface domain's cell volume.
    - `"molar density`" Defaults to domain's molar_density_liquid.




Top cell
^^^^^^^^


Arbitrary function
^^^^^^^^^^^^^^^^^^
 A secondary variable evaluator which evaluates functions on its dependenecies.
Uses functions to evaluate arbitrary algebraic functions of its dependencies.

For example, one might write a dependency:

  VARNAME = 0.2 * DEP1 - DEP2 + 3

as:

Note this is not done by region currently, but could easily be extended to do
so if it was found useful.

Example:
..xml:
    <ParameterList name="VARNAME">
      <Parameter name="field evaluator type" type="string" value="secondary variable from function"/>
      <Parameter name="evaluator dependencies" type="Array{string}" value="{DEP1, DEP2}"/>
      <ParameterList name="function">
        <ParameterList name="function-linear">
          <Parameter name="x0" type="Array(double)" value="{0.0,0.0}" />
          <Parameter name="y0" type="double" value="3." />
          <Parameter name="gradient" type="Array(double)" value="{0.2, -1}" />
        </ParameterList>
      </ParameterList>
    </ParameterList>

 



Not-so-generic evaluators
-------------------------
Nearly all of these probably can or should be refactored into a more
generic evaluator.  If not, they should get moved to a more-correct
section of the documentation.

Column Average Temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^
Can this become an option in `Column Summation`_?




Column Water Content
^^^^^^^^^^^^^^^^^^^^
Can this become a `Column Summation`_?



Gas Content
^^^^^^^^^^^
Can this become a `Column Summation`_?



Max Thaw Depth
^^^^^^^^^^^^^^
Can this become a generic "Maximum" or similar?



Thaw Depth (Columns)
^^^^^^^^^^^^^^^^^^^^
Can this become a `Subgrid Aggregation`_ of a `Thaw Depth`_?



Water Table (Columns)
^^^^^^^^^^^^^^^^^^^^^
Can this become a `Subgrid Aggregation`_  of a `Water Table`_?






    

InitialConditions
=================

Initial condition specs are used in three places:

* In the `"initial conditions`" sublist of state, in which the value
  of atomic constants are provided (not really initial conditions and
  should be renamed).  These atomic values are not controlled by
  evaluators, and are not included in the DaG.

* Within the PK_ spec which describes the initial condition of primary variables (true
  initial conditions).

* In `IndependentVariableEvaluator Constant <Constant_>`_

The first may be of multiple types of data, while the latter two are
nearly always fields on a mesh (e.g. CompositeVectors).  The specific
available options differ by data type.

Scalars
-------


``[initial-conditions-scalar-spec]``

* `"value`" ``[double]`` Set the value.


Example:

.. code-block:: xml

  <ParameterList name="atmospheric_pressure">
    <Parameter name="value" type="double" value="101325."/>
  </ParameterList>





Vectors of Scalars
------------------


``[initial-conditions-constantvector-spec]``

* `"value`" ``[Array(double)]`` Set the value.


Example:

.. code-block:: xml

  <ParameterList name="gravity">
    <Parameter name="value" type="Array(double)" value="{{0.0, -9.81}}"/>
  </ParameterList>




Fields on a Mesh
----------------


The vast majority of Fields are vectors of data on a mesh entity.  This is a
CompositeVector in Amanzi-ATS (as it may consist of multiple components, each
on its own mesh entity).  Initializing these components may be done in avariety
of ways:


``[initial-conditions-spec]``

* `"restart file`" ``[string]`` **optional** If provided, read IC value from a
  checkpoint file of this name.

* `"cells from file`" ``[string]`` **optional** If provided, read IC value for
  cells from a checkpoint file of this name.

* `"exodus file initialization`" ``[exodus-file-reader-spec]`` **optional** If provided, initialize data from an Exodus file.


* `"constant`" ``[double]`` **optional** If provided, set the IC to this value everywhere.

* `"initialize from 1D column`" ``[initialize-1d-column-spec]`` **optional**
   If provided, initialize data by interpolating from a 1D solution in depth
   relative to the surface.

* `"function`" ``[composite-vector-function-spec-list]`` **optional** If
   provided, a region-based list of functions to evaluate piecewise over the
   domain.

* `"initialize faces from cell`" ``[bool]`` **false** In the case of mimetic
   or other face-inclusive discretizations, calculate initial face data by
   averages of the neighboring cells.

* `"write checkpoint`" ``[bool]`` **true** Write this data when checkpointing.

* `"write vis`" ``[bool]`` **true** Write this data when visualizing.


Capability for 1D column initialization:

``[initialize-1d-column-spec]``
* `"file`" ``[string]`` Filename including data.
* `"z header`" ``[string]`` **"/z"**  Name of the data in the above file for z-coordinate of this field.
* `"f header`" ``[string]`` **"/field-name"** Name of the data in the above file for this field.
* `"coordinate orientation`" ``[string]`` **"standard"** Either `"standard`" or `"depth`" -- is positive z up (standard) or down (depth)?

ONE OF:
* `"surface sideset`" ``[string]`` Name of the face-set in the mesh which determines the zero surface for interpolation.  This must cover the entire top surface boundary.
OR:
* `"surface sidesets`" ``[Array(string)]`` List of the names for the case of multiple sidesets needed.
END

ExodusII file reader control:

``[exodus-file-reader-spec]``
* `"file`" ``[string]`` Filename
* `"attributes`" ``[Array(string)]`` List of Exodus attributes to read, unclear when/if/how this works?





BoundaryConditions
===================



In general, boundary conditions are provided in a hierarchical list by
boundary condition type, then functional form.  Boundary condition specs are
split between two types -- those which require a user-provided function
(i.e. Dirichlet data, etc) and those which do not (i.e. zero gradient
conditions).

A list of conditions might pull in both Dirichlet and Neumann data on
different regions, or use different functions on different regions.  The
following example illustrates how boundary conditions are prescribed across
the domain for a typical PK:

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="DIRICHLET_TYPE">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
     <ParameterList name="BC east">
       <Parameter name="regions" type="Array(string)" value="{east}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="102325."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="mass flux">
     <ParameterList name="BC north">
       <Parameter name="regions" type="Array(string)" value="{north}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="zero gradient">
     <ParameterList name="BC south">
       <Parameter name="regions" type="Array(string)" value="{south}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Different PKs populate this general format with different names, replacing
DIRICHLET_TYPE and DIRICHLET_FUNCTION_NAME.
  
 


Flow-specific Boundary Conditions
----------------------------------



Flow boundary conditions must follow the general format shown in
BoundaryConditions_.  Specific conditions implemented include:

Dirichlet (pressure) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for both surface and subsurface flows, this provides pressure data on
boundaries (in [Pa]).

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="pressure">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary pressure">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Dirichlet (head) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Used for both surface and subsurface flows, this provides head data (in [m]
above the land surface), typically as a function of x & y.  In the subsurface
case, the z-value is given by hydrostatic relative to that head.

.. math::
  p = p_{atm} + rho * g * (h(x,y) + z_{surf} - z)

where h is the head function provided.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="head">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary head">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.01"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Dirichlet (fixed level) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This fixes the water table at a constant elevation.  It is a head condition
that adapts to the surface elevation, adjusting the head to a datum that is a
fixed absolute z coordinate.

.. math::
  p = p_{atm} + rho * g * (h(x,y) - z)

where h is the head function provided.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="fixed level">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="fixed level">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Neumann (mass flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for both surface and subsurface flows, this provides mass flux data (in [mol m^-2 s^-1], for the subsurface, or [mol m^-1 s^-1] for the surface, in the outward normal direction) on boundaries.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="mass flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="-1.e-3"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Neumann (fix level flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for surface only,  this provides fixed level ([m])  velocity data (in [mol m^-3 s^-1], in the outward normal direction on boundaries.

Example:

.. code-block:: xml
    
     <ParameterList name="boundary conditions">
       <ParameterList name="fixed level flux">
          <ParameterList name="river level south">
            <Parameter name="regions" type="Array(string)" value="{river south}"/>
            <ParameterList name="fixed level">
               <ParameterList name="function-constant">
                 <Parameter name="value" type="double" value="0.5"/>
               </ParameterList>
            </ParameterList>
            <ParameterList name="velocity">
               <ParameterList name="function-constant">
                 <Parameter name="value" type="double" value="-25"/>
               </ParameterList>
            </ParameterList>            
          </ParameterList>
       </ParameterList>
    </ParameterList>


 
Seepage face boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A variety of seepage face boundary conditions are permitted for both surface
and subsurface flow PKs.  Typically seepage conditions are of the form:

  - if :math:`q \cdot \hat{n} < 0`, then :math:`q = 0`
  - if :math:`p > p0`, then :math:`p = p0`

This ensures that flow is only out of the domain, but that the max pressure on
the boundary is specified by :math:`p0`.

Example: pressure (for surface or subsurface)

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face pressure">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary pressure">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Example: head (for surface)
 
.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face head">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary head">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Additionally, an infiltration flux may be prescribed, which describes the max
flux.  This is for surface faces on which a typical precipitation rate might
be prescribed, to be enforced until the water table rises to the surface, at
which point the precip is turned off and water seeps into runoff.  This
capability is experimental and has not been well tested.

  - if :math:`q \cdot \hat{n} < q_0`, then :math:`q = q_0`
  - if :math:`p > p_{atm}`, then :math:`p = p_{atm}`

Note the condition also accepts a parameter:

* `"explicit time index`" ``[bool]`` **false** If true, the _type_ of the BC is
  evaluated at the old time, keeping it fixed while the nonlinear solve
  iterates.

Example: seepage with infiltration

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face with infiltration">
     <ParameterList name="BC west">
       <Parameter name="explicit time index" type="bool" value="true"/>
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="-1.e-5"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Note it would be straightforward to add both p0 and q0 in the same condition;
this has simply not had a use case yet.


Zero head gradient boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for surface flows, this is an "outlet" boundary condition which looks to
enforce the condition that

.. math::
  \div h \cdot \hat{n} = 0

for head :math:`h` and outward normal :math:`\hat{n}`.  Note that this is an
"outlet" boundary, in the sense that it should really not be used on a
boundary in which

.. math::
  \div z \cdot \hat{n} > 0.

This makes it a useful boundary condition for benchmark and 2D problems, where
the elevation gradient is clear, but not so useful for DEM-based meshes.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="zero gradient">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Critical depth boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Also for surface flows, this is an "outlet" boundary condition which looks to
set an outward flux to take away runoff.  This condition is given by:

.. math::
  q = \sqrt{g \hat{z}} n_{liq} h^1.5

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="critical depth">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Dynamic boundary condutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The type of boundary conditions maybe changed in time depending on the switch function of TIME.

.. code-block:: xml

   <ParameterList name="dynamic">
     <Parameter name="regions" type="Array(string)" value="{surface west}"/>
     <ParameterList name="switch function">
       <ParameterList name="function-tabular">
         <Parameter name="file" type="string" value="../data/floodplain2.h5" />
         <Parameter name="x header" type="string" value="Time" />
         <Parameter name="y header" type="string" value="Switch" />
         <Parameter name="form" type="Array(string)" value="{constant}"/>
       </ParameterList>
     </ParameterList>
          
     <ParameterList name="bcs">
       <Parameter name="bc types" type="Array(string)" value="{head, mass flux}"/>
       <Parameter name="bc functions" type="Array(string)" value="{boundary head, outward mass flux}"/>

       <ParameterList name="mass flux">
         <ParameterList name="BC west">
           <Parameter name="regions" type="Array(string)" value="{surface west}"/>
           <ParameterList name="outward mass flux">
             <ParameterList name="function-tabular">
               <Parameter name="file" type="string" value="../data/floodplain2.h5" />
               <Parameter name="x header" type="string" value="Time" />
               <Parameter name="y header" type="string" value="Flux" />
               <Parameter name="form" type="Array(string)" value="{linear}"/>
             </ParameterList>
            </ParameterList>
          </ParameterList>
       </ParameterList>

       <ParameterList name="head">  
          <ParameterList name="BC west">
            <Parameter name="regions" type="Array(string)" value="{surface west}"/>
            <ParameterList name="boundary head">
              <ParameterList name="function-tabular">
                 <Parameter name="file" type="string" value="../data/floodplain2.h5" />
                 <Parameter name="x header" type="string" value="Time" />
                 <Parameter name="y header" type="string" value="Head" />
                 <Parameter name="form" type="Array(string)" value="{linear}"/>
               </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
     </ParameterList>

   </ParameterList>

 


Transport-specific Boundary Conditions
--------------------------------------




Energy-specific Boundary Conditions
-----------------------------------








Time integrators, solvers, and other mathematical specs
#######################################################

  Common specs for all solvers and time integrators, used in PKs.

There are three commonly used broad classes of time integration
strategies.

"Update" methods are the simplest -- they use no formal
mathematical definition of a differential equation, but instead
implicitly use a process by which variables at the new time are
directly calculated.  Typically there *is* an implied ODE or PDE here,
but it is not stated as such and time integration routines are not
used.  Examples of these are common in biogeochemistry and vegetation
models.

"Explicit" time methods are the next simplest.  These include a
variety of options from forward Euler to higher order Runge-Kutta
schemes.  These only require evaluating forward models where we have
existing of the dependencies.  If they work, these are great thanks to
their deterministic nature and lack of expensive, memory-bandwith
limited solvers.  But they only work on some types of problems.
Examples of of these include transport, where we use high order time
integration schemes to preserve fronts.

"Implicit" and semi-implicit methods instead require the evaluation of
a residual equation -- the solution is guessed at, and the residual is
calculated, which measures how far the equation is from being
satisfied.  This measure is then inverted, finding a correction to the
guess which hopefully reduces the residual.  As the residual goes to
zero, the error, a measure of the difference between the guess and the
true solution, also goes to zero.  To do this inversion, we lean on
Newton and Newton-like methods, which attempt to somehow linearize, or
approximately linearize, the residual function near the guess in order
to calculate an update.  In this case, the time integration scheme
requires both a nonlinear solver (to drive the residual to zero) and a
linear solver or approximate solver (to calculate the correction).

TimeIntegrator
==============

Currently there are two classes of time integration schemes used in
ATS: explicit (including a range of single and multi-stage) methods
and BDF1, or Backward Euler.

Explicit Time Integration
-------------------------
 Explicit time integration methods in a generalized form.

This class implements several explicit Runge Kutta methods:

- forward Euler     (1st order)  --> "forward_euler"
- Heun-Euler method (2nd order)  --> "heun_euler"
- Midpoint method   (2nd order)  --> "midpoint"
- Ralston method    (2nd order)  --> "ralston"
- TVD RK method     (3rd order)  --> "tvd_3rd_order"
- Kutta method      (3rd order)  --> "kutta_3rd_order"
- Runge Kutta       (4th order)  --> "runge_kutta_4th_order"
- User defined      (whatever)   --> user_defined, use special constructor to create

Note that user-defined is only for developers currently, and cannot be created from an input file.

The RK tableau is made up of the three private objects a, b, and c below.  they
are arranged as follows:

.. code-block:

    c[0]   | 
    c[1]   | a(1,0)
     .     | a(2,0)    a(2,1)
     .     |   .         .        
     .     |   .         . 
    c[s-1 ]| a(s-1,0)  a(s-1,1)  . . .  a(s-1,s-2) 
    ---------------------------------------------------------
           |   b[0]      b[1]    . . .    b[s-2]      b[s-1] 
  
Note that c[0] should always equal zero, and that the entries in the matrix a
that are not listed in this tableau are not used
  
The implemented general Runge Kutta scheme of order s based on this tableau arrangement is

.. math::
    y_{n+1} = y_n + \sum{i=0}^{s-1} b[i]*k_i
  
    with 
  
      k_0 = h * f(t_n, y_n) \\
      k_1 = h * f(t_n + c[1]*h, y_n + a(1,0)*k_0) \\
      k_2 = h * f(t_n + c[2]*h, y_n + a(2,0)*k_0 + a(2,1)*k_1) \\
       . \\
       . \\
       . \\ 
      k_{s-1} = h * f(t_n + c[s-1]*h, y_n + a(s-1,0)*k_0 + ... + a(s-1,s-2)*k_{s-2})


.. _explicit-ti-rk-spec:
.. admonition:: explicit-ti-rk-spec

    * `"verbose object`" ``[verbose-object-spec]`` A `Verbose Object`_

    * `"RK method`" ``[string]`` **forward euler**  One of: `"forward Euler`", `"heun euler`", `"midpoint`", `"ralston`", `"tvd 3rd order`", `"kutta 3rd order`", `"runge kutta 4th order`"
      



Backward Euler
--------------
 Solves globally implicit systems using backward Euler

Backward Euler is the simplest of the implicit methods.  It solves time
integration schemes by evaluating all time derivatives at the new time.  This
makes it unconditionally stable, though potentially not very accurate.  This
unconditional stability tends to make it the workhorse of the types of stiff,
nonlinear parabolic equations such as Richards equation and the diffusion wave
approximation.

In this method, we look to solve:

.. math::
    \frac{\partial \mathbf{u}}{\partial t} = f(\mathbf{u},\mathbf{x},t)

via the time discretization scheme:

.. math::
    \frac{\mathbf{u}^{t + \Delta t} - \mathbf{u}^{t}}{\Delta t} = f(\mathbf{u}^{t + \Delta t}, \mathbf{x}, t + \Delta t)

.. _bdf1-ti-spec:
.. admonition:: bdf1-ti-spec

    * `"verbose object`" ``[verbose-object-spec]`` A `Verbose Object`_

    * `"residual debugger`" ``[residual-debugger-spec]`` A `Residual Debugger`_ object.

    * `"max preconditioner lag iterations`" ``[int]`` **0** specifies frequency
      of preconditioner recalculation.

    * `"freeze preconditioner`" ``[bool]`` **false** enforces preconditioner to
      be updated only once per non-linear solver. When set to true, the above
      parameter is ignored.

    * `"extrapolate initial guess`" ``[bool]`` **true** identifies forward time
      extrapolation of the initial guess.

    * `"nonlinear iteration initial guess extrapolation order`" ``[int]`` **1**
      defines extrapolation algorithm. Zero value implies no extrapolation.

    * `"restart tolerance relaxation factor`" ``[double]`` **1** Changes the
      nonlinear tolerance on restart. The time integrator is usually restarted
      when a boundary condition changes drastically. It may be beneficial to
      loosen the nonlinear tolerance on the first several time steps after the
      time integrator restart. The default value is 1, while a reasonable value
      may be as large as 1000.

    * `"restart tolerance relaxation factor damping`" ``[double]`` **1**
      Controls how fast the loosened nonlinear tolerance will revert back to
      the one specified in `"nonlinear tolerance`". If the nonlinear tolerance
      is "tol", the relaxation factor is "factor", and the damping is "d", and
      the time step count is "n" then the actual nonlinear tolerance is "tol *
      max(1.0, factor * d ** n)". Reasonable values are between 0 and 1.

    INCLUDES
    - ``[solver-typed-spec]`` *Uses a* Solver_.
    - ``[timestep-controller-typed-spec]`` *Uses a* `Timestep Controller`_


Note this also accepts an object that provides the `BDF1 Solver Interface`_.




BDF1 Solver Interface
^^^^^^^^^^^^^^^^^^^^^


Timestep Controller
-------------------
 Factory for creating TimestepController objects

A TimestepController object sets what size timestep to take.  This can be a
variety of things, from fixed timestep size, to adaptive based upon error
control, to adapter based upon simple nonlinear iteration counts.

Available types include:

- `Timestep Controller Fixed`_  (type `"fixed`"), a constant timestep
- `Timestep Controller Standard`_ (type `'standard`"), an adaptive timestep based upon nonlinear iterations
- `Timestep Controller Smarter`_ (type `'smarter`"), an adaptive timestep based upon nonlinear iterations with more control
- `Timestep Controller Adaptive`_ (type `"adaptive`"), an adaptive timestep based upon error control.
- `Timestep Controller From File`_ (type `"from file`"), uses a timestep history loaded from an HDF5 file.  (Usually only used for regression testing.)


.. _timestep-controller-typed-spec:
.. admonition:: timestep-controller-typed-spec

    * `"timestep controller type`" ``[string]`` Set the type.  One of: `"fixed`", `"standard`", `"smarter`", `"adaptive`", or `"from file`"
    * `"timestep controller X parameters`" ``[list]`` List of parameters for a timestep controller of type X.




Timestep Controller Fixed
^^^^^^^^^^^^^^^^^^^^^^^^^
  Timestep controller providing constant timestep size.

``TimestepControllerFixed`` is a simple timestep control mechanism which sets
a constant timestep size.  Note that the actual timestep size is given by the
minimum of PK's initial timestep sizes.

No parameters are required.




Timestep Controller Standard
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Simple timestep control based upon previous iteration count.

This is a simple timestep control mechanism
which sets the next timestep based upon the previous timestep and how many
nonlinear iterations the previous timestep took to converge.

The timestep for step :math:`k+1`, :math:`\Delta t_{k+1}`, is given by:

- if :math:`N_k > N^{max}` then :math:`\Delta t_{k+1} = f_{reduction} * \Delta t_{k}`
- if :math:`N_k < N^{min}` then :math:`\Delta t_{k+1} = f_{increase} * \Delta t_{k}`
- otherwise :math:`\Delta t_{k+1} = \Delta t_{k}`

where :math:`\Delta t_{k}` is the previous timestep and :math:`N_k` is the number of 
nonlinear iterations required to solve step :math:`k`:.

.. _timestep-controller-standard-spec:
.. admonition:: timestep-controller-standard-spec

    * `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the previous step took more than this.
    * `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the previous step took less than this.
    * `"time step reduction factor`" ``[double]`` :math:`f_{reduction}`, reduce the previous timestep by this multiple.
    * `"time step increase factor`" ``[double]`` :math:`f_{increase}`, increase the previous timestep by this multiple.
    * `"max time step`" ``[double]`` The max timestep size allowed.
    * `"min time step`" ``[double]`` The min timestep size allowed.  If the step has failed and the new step is below this cutoff, the simulation fails.




Timestep Controller Smarter
^^^^^^^^^^^^^^^^^^^^^^^^^^^
  Slightly smarter timestep controller based upon a history of previous timesteps.

This is based on `Timestep Controller Standard`_, but also tries to be a bit
smarter to avoid repeated increase/decrease loops where the step size
decreases, converges in few iterations, increases, but then fails again.  It
also tries to grow the step geometrically to more quickly recover from tricky
nonlinearities.

.. _timestep-controller-smarter-spec:
.. admonition:: timestep-controller-smarter-spec

    * `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the previous step took more than this.
    * `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the previous step took less than this.
    * `"time step reduction factor`" ``[double]`` :math:`f_{reduction}`, reduce the previous timestep by this multiple.
    * `"time step increase factor`" ``[double]`` :math:`f_{increase}`, increase the previous timestep by this multiple.  Note that this can be modified geometrically in the case of repeated successful steps.
    * `"max time step increase factor`" ``[double]`` **10.** The max :math:`f_{increase}` will ever get.
    * `"growth wait after fail`" ``[int]`` Wait at least this many timesteps before attempting to grow the timestep after a failed timestep.
    * `"count before increasing increase factor`" ``[int]`` Require this many successive increasions before multiplying :math:`f_{increase}` by itself.





Timestep Controller Adaptive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Adaptive timestep control based upon previous iteration count.

This is under development and is based on a posteriori error estimates.




Timestep Controller From File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  Timestep controller which loads a timestep history from file.

This loads a timestep history from a file, then advances the step size with
those values.  This is mostly used for testing purposes, where we need to force
the same timestep history as previous runs to do regression testing.  Otherwise
even machine roundoff can eventually alter number of iterations enough to alter
the timestep history, resulting in solutions which are enough different to
cause doubt over their correctness.

.. _timestep-controller-from-file-spec:
.. admonition:: timestep-controller-from-file-spec

    * `"file name`" ``[string]`` Path to hdf5 file containing timestep information.
    * `"timestep header`" ``[string]`` Name of the dataset containing the history of timestep sizes.





Nonlinear Solver
================
.. _Solver:
   
 A factory for creating nonlinear solvers.

Nonlinear solvers are used within implicit time integration schemes to drive
the residual to zero and thereby solve for the primary variable at the new
time.

.. _solver-typed-spec:
.. admonition:: solver-typed-spec

    * `"solver type`" ``[string]`` Type of the solver.  One of:

      - `"Newton`" See `Solver: Newton and Inexact Newton`_
      - `"JFNK`" See `Solver: Jacobian-Free Newton Krylov`_
      - `"line search`" See `Solver: Newton with Line Search`_
      - `"continuation`" See `Solver: Nonlinear Continuation`_
      - `"nka`" See `Solver: Nonlinear Krylov Acceleration`_
      - `"aa`" See `Solver: Anderson Acceleration`_
      - `"nka line search`" See `Solver: NKA with Line Search`"
      - `"nka_ls_ats`" See `Solver: NKA with Line Search, ATS`_
      - `"nka_bt_ats`" See `Solver: NKA with backtracking, ATS`_
      - `"nox`" See `Solver: NOX`_
    
    * `"_solver_type_ parameters`" ``[_solver_type_-spec]`` A sublist containing
      parameters specific to the type.

.. warning::

    `"JFNK`", `"line search`", and `"continuation`" methods have not been
    beaten on as much as other methods.  `"nka_ls_ats`" is somewhat deprecated
    and probably shouldn't be used.  Prefer `"nka`" for simple problems,
    `"nka_bt_ats`" for freeze-thaw problems or other problems with strong
    nonlinearities, and `"Newton`" when you have a good Jacobian.  While
    `"nox`" hasn't been used extensively, it may be quite useful.


      



Solver: Newton and Inexact Newton
---------------------------------
 Straightforward Newton/Inexact Newton solver.

The classical Newton method works well for cases where Jacobian is available
and corresponds to a stable (e.g. upwind) discretization.  The inexact Newton
methods work for cases where the discrete Jacobian is either not available, or
not stable, or computationally expensive. The discrete Jacobian is replaced by
a stable approximation of the continuum Jacobian. The choice between exact and
inexact is not made by the Solver, but instead by the PK.  Both use the
ApplyPreconditioner() method -- if this applies the true Jacobian, then the
method is Newton.  If it applies an appoximation, it is inexact Newton.

.. _solver-newton-spec:
.. admonition:: solver-newton-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** defines the required error
      tolerance. The error is calculated by a PK.
    
    * `"monitor`" ``[string]`` **monitor update** specifies control of the
      nonlinear residual. The available options are `"monitor update`" and
      `"monitor residual`".

    * `"limit iterations`" ``[int]`` **50** defines the maximum allowed number
      of iterations.

    * `"diverged tolerance`" ``[double]`` **1.e10** defines the error level
      indicating divergence of the solver. The error is calculated by a PK.

    * `"max du growth factor`" ``[double]`` **1.e5** allows the solver to
      identify divergence pattern on earlier iterations. If the maximum norm of
      the solution increment changes drastically on two consecutive iterations,
      the solver is terminated.

    * `"max error growth factor`" ``[double]`` **1.e5** defines another way to
      identify divergence pattern on earlier iterations. If the PK-specific
      error changes drastically on two consecutive iterations, the solver is
      terminated.

    * `"max divergent iterations`" ``[int]`` **3** defines another way to
      identify divergence pattern on earlier iterations. If the maximum norm of
      the solution increment grows on too many consecutive iterations, the
      solver is terminated.

    * `"modify correction`" ``[bool]`` **true** allows a PK to modify the
      solution increment. One example is a physics-based clipping of extreme
      solution values.

    * `"stagnation iteration check`" ``[int]`` **8** determines the number of
      iterations before the stagnation check is turned on. The stagnation
      happens when the current L2-error exceeds the initial L2-error.

 


Solver: Jacobian-Free Newton Krylov
-----------------------------------
 Decorator for using a Solver with JFNK as the preconditioner.

Jacobian-Free Newton Krylov uses a finite difference scheme to approximate the
action of the Jacobian matrix, then uses a Krylov method (which only needs the
action of the Jacobian and not the Jacobian itself) to calculate the action of
the inverse of the Jacobian, thereby providing a Newton-like update.  As the
linear Krylov scheme converges to the inverse action, the nonlinear solution
converges to the same solution as a true Newton method.

This implementation simply replaces a SolverFnBase's ApplyPreconditioner() with
a new ApplyPreconditioner() which uses the Krylov method with the action of the
forward operator to (hopefully) improve, relative to the supplied approximate
inverse, the estimate of the inverse.

.. _solver-jfnk-spec:
.. admonition:: solver-jfnk-spec

    * `"nonlinear solver`" ``[solver-typed-spec]`` The outer nonlinear solver to use.

    * `"inverse`" ``[inverse-typed-spec]`` The Krylov method to use.

    * `"JF matrix parameters`" ``[jf-matrix-spec]`` See jf-matrix-spec_





The Jacobian-Free Matrix operator, which is used to estimate the action of the
Jacobian.

A variety of methods are available for choosing the epsilon used to approximate
the action of the Jacobian.  They are documented in Knoll & Keyes 2004 paper.


..todo:: Document these

.. _jf-matrix-spec:
.. admonition:: jf-matrix-spec

    * `"typical solution value`" ``[double]`` **100** Used in relative action
      approximations. OPTION NOT IMPLEMENTED

    * `"finite difference epsilon`" ``[double]`` **1.e-8** defines the base finite
      difference epsilon.

    * `"method for epsilon`" ``[string]`` defines a method for calculating finite
      difference epsilon. Available option is "Knoll-Keyes", "Knoll-Keyes L2",
      "Brown-Saad".  See Knoll




Solver: Newton with Line Search
-------------------------------
 Backtracking line search on the provided correction as a solver.

Line Search accepts a correction from the Jacobian, then uses a
process to attempt to minimize or at least ensure a reduction in the residual
while searching *in that direction*, but not necessarily with the same
magnitude, as the provided correction.  The scalar multiple of the search
direction is given by :math:`\alpha`.

This globalization recognizes that a true inverse Jacobian is a local
measurement of the steepest descent direction, and so while the direction is
guaranteed to be the direction which best reduces the residual, it may not
provide the correct magnitude.

Note, this always monitors the residual.

.. _solver-backtracking-spec:
.. admonition:: solver-backtracking-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** defines the required error
      tolerance. The error is calculated by a PK.
    
    * `"limit iterations`" ``[int]`` **50** defines the maximum allowed number
      of iterations.

    * `"diverged tolerance`" ``[double]`` **1.e10** defines the error level
      indicating divergence of the solver. The error is calculated by a PK.
      Set to a negative value to ignore this check.

    * `"max error growth factor`" ``[double]`` **1.e5** defines another way to
      identify divergence pattern on earlier iterations. If the PK-specific
      error changes drastically on two consecutive iterations, the solver is
      terminated.

    * `"modify correction`" ``[bool]`` **false** allows a PK to modify the
      solution increment. One example is a physics-based clipping of extreme
      solution values.

    * `"accuracy of line search minimum [bits]`" ``[int]`` **10**

    * `"min valid alpha`" ``[double]`` **0** 

    * `"max valid alpha`" ``[double]`` **10.**

    * `"max line search iterations`" ``[int]`` **10** 

  
 


Solver: Nonlinear Continuation
------------------------------
 A very simple nonlinear continuation method.

Continuation methods are useful when the nonlinearity can be controlled by a
single simple parameter.  In this method, the nonlinear problem is solved with
a less-nonlinear value of the parameter, and the solution of that is used as
the initial guess to solve a harder problem.  As each successive problem is
solved, the continuation parameter is changed closer and closer to the true
value.

Few if any PKs support this method currently -- it requires the PK to provide more
interface about how to update the continuation parameter.

.. _solver-continuation-spec:
.. admonition:: solver-continuation-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** defines the required error
      tolerance. The error is calculated by a PK.
    
    * `"number of continuation steps`" ``[int]`` **5** How many steps to take
      from initial parameter to final parameter.

    * `"inner solver`" ``[solver-typed-spec]``A Solver_, used at each step.

 


Solver: Nonlinear Krylov Acceleration
-------------------------------------
 Nonlinear Krylov Acceleration as a nonlinear solver.

Uses the Nonlinear Krylov acceleration method of Carlson and Miller to do
effectively a multivariant secant method, accelerating the solution of a
nonlinear solve.  This method can be significantly faster than Newton,
especially with an approximate Jacobian.

  Calef et al. "Nonlinear Krylov acceleration applied to a discrete ordinates
  formulation of the k-eigenvalue problem." JCP 238 (2013): 188-209.

  N. N. Carlson, K. Miller, Design and application of a gradient-weighted
  moving finite element code II: In two dimensions, SIAM J. Sci.  Comput. 19
  (3) (1998) 766–798.


.. _solver-nka-spec:
.. admonition:: solver-nka-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** Defines the required error
      tolerance. The error is calculated by a PK.

    * `"monitor`" ``[string]`` **monitor update** Specifies control of the
      nonlinear residual. The available options are `"monitor update`",
      `"monitor residual`", `"monitor preconditioned residual`", `"monitor l2
      residual`", and `"monitor preconditioned l2 residual`".

    * `"limit iterations`" ``[int]`` **20** Defines the maximum allowed number
      of iterations.

    * `"diverged tolerance`" ``[double]`` **1.e10** Defines the error level
      indicating divergence of the solver. The error is calculated by a PK.
      Set to a negative value to ignore this check.

    * `"diverged l2 tolerance`" ``[double]`` **1.e10** Defines another way to
      identify divergence of the solver. If the relative l2 (little l) norm of
      the solution increment is above this value, the solver is terminated.
      Set to a negative value to ignore this check.

    * `"diverged pc tolerance`" ``[double]`` **1e10** Defines another way to
      identify divergence of the solver. If the relative maximum norm of the
      solution increment (with respect to the initial increment) is above this
      value, the solver is terminated.
      Set to a negative value to ignore this check.

    * `"diverged residual tolerance`" ``[double]`` **1e10** Defines another way
      to identify divergence of the solver. If the relative l2 norm of the
      residual (with respect to the initial residual) is above this value, the
      solver is terminated.  Set to a negative value to ignore this check.

    * `"max du growth factor`" ``[double]`` **1e5** Allows the solver to
      identify divergence pattern on earlier iterations. If the maximum norm of
      the solution increment changes drastically on two consecutive iterations,
      the solver is terminated.

    * `"max error growth factor`" ``[double]`` **1e5** Defines another way to
      identify divergence pattern on earlier iterations. If the PK-specific
      error changes drastically on two consecutive iterations, the solver is
      terminated.

    * `"max divergent iterations`" ``[int]`` **3** Defines another way to
      identify divergence pattern on earlier iterations. If the maximum norm of
      the solution increment grows on too many consecutive iterations, the
      solver is terminated.

    * `"modify correction`" ``[bool]`` **false** Allows a PK to modify the
      solution increment. One example is a physics-based clipping of extreme
      solution values.

    * `"lag iterations`" ``[int]`` **0** Delays the NKA acceleration, but
      updates the Krylov space.

    * `"max nka vectors`" ``[int]`` **10** Defines the maximum number of
      consecutive vectors used for a local space.

    * `"nka vector tolerance`" ``[double]`` **0.05** Defines the minimum
      allowed orthogonality between vectors in the local space. If a new vector
      does not satisfy this requirement, the space is modified.

 


Solver: Anderson Acceleration
-----------------------------
 Anderson acceleration as a nonlinear solver.

This is a variation of the GMRES solver for nonlinear problems.

.. _solver-aa-spec:
.. admonition:: solver-aa-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** Defines the required error
      tolerance. The error is calculated by a PK.

    * `"limit iterations`" ``[int]`` **20** Defines the maximum allowed number
      of iterations.

    * `"diverged tolerance`" ``[double]`` **1.e10** Defines the error level
      indicating divergence of the solver. The error is calculated by a PK.
      Set to a negative value to ignore this check.

    * `"diverged l2 tolerance`" ``[double]`` **1.e10** Defines another way to
      identify divergence of the solver. If the relative L2 norm of the
      solution increment is above this value, the solver is terminated.
      Set to a negative value to ignore this check.

    * `"max du growth factor`" ``[double]`` **1e5** Allows the solver to
      identify divergence pattern on earlier iterations. If the maximum norm of
      the solution increment changes drastically on two consecutive iterations,
      the solver is terminated.

    * `"max divergent iterations`" ``[int]`` **3** Defines another way to
      identify divergence pattern on earlier iterations. If the maximum norm of
      the solution increment grows on too many consecutive iterations, the
      solver is terminated.
      
    * `"max aa vectors`" ``[int]`` **10** Defines the maximum number of
      consecutive vectors used for a local space.

    * `"modify correction`" ``[bool]`` **false** Allows a PK to modify the
      solution increment. One example is a physics-based clipping of extreme
      solution values.

    * `"relaxation parameter`" ``[double]`` **1** Damping factor for increment.

 


Solver: NKA with Line Search
----------------------------
 NKA nonlinear solver with a line-search based on a Brendt minimization algorithm.

Does NKA, then checks if that correction has reduced the residual by at least a
tolerance.  If not, this uses a Brendt minimization algorithm to try and find
an :math:`\alpha` that minimizes the reduction in the residual.

Note, this always monitors the residual.

.. _solver-nka-ls-spec:
.. admonition solver-nka-ls-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** Defines the required error
      tolerance. The error is calculated by a PK.

    * `"limit iterations`" ``[int]`` **20** Defines the maximum allowed number
      of iterations.

    * `"backtrack monitor`" ``[string]`` **monitor either** What norm is
      checked to determine whether backtracking has improved the residual or
      not?  One of `"monitor enorm`", `"monitor L2 residual`", or `'monitor
      either`"

    * `"backtrack tolerance`" ``[double]`` **0.** If the default update reduces
      the residual by at least this much, line search is not performed.

    * `"accuracy of line search minimum [bits]`" ``[int]`` **10** Convergence criteria on Brendt algorithm.

    * `"min valid alpha`" ``[double]`` **0** Lower bound on Brendt algorithm.

    * `"max valid alpha`" ``[double]`` **10.** Upper bound on Brendt algorithm.

    * `"max line search iterations`" ``[int]`` **10** Max iterations for the Brendt algorithm.

    * `"max nka vectors`" ``[int]`` **10** Defines the maximum number of
      consecutive vectors used for a local space.

    * `"nka vector tolerance`" ``[double]`` **0.05** Defines the minimum
      allowed orthogonality between vectors in the local space. If a new vector
      does not satisfy this requirement, the space is modified.

  



Solver: NKA with Line Search, ATS
---------------------------------


Solver: NKA with backtracking, ATS
----------------------------------
 Nonlinear solve using NKA with a heuristic based backtracking.

Whereas line search uses a formal minimization method, backtracking simply uses
a heuristic multiplier on :math:`\alpha` to find a correction that sufficiently
reduces the residual.  This can be significantly faster than the full
minimization problem, and finding the true minimum may not be as important as
simply doing better and going on to the next nonlinear iteration.

This is the workhorse for hard ATS problems, as it is usually rather efficient,
even in problems where the linear solve results in a correction that is way too
large (e.g. for steep nonlinearities such as phase change).

Note this always monitors the residual, and the correction is always modified.

.. _solver-nka-bt-ats-spec:
.. admonition:: solver-nka-bt-ats-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** Defines the required error
      tolerance. The error is calculated by a PK.

    * `"limit iterations`" ``[int]`` **20** Defines the maximum allowed number
      of iterations.

    * `"diverged tolerance`" ``[double]`` **1.e10** Defines the error level
      indicating divergence of the solver. The error is calculated by a PK.
      Set to a negative value to ignore this check.

    * `"nka lag iterations`" ``[int]`` **0** Delays the NKA acceleration, but
      updates the Krylov space.

    * `"max nka vectors`" ``[int]`` **10** Defines the maximum number of
      consecutive vectors used for a local space.

    * `"nka vector tolerance`" ``[double]`` **0.05** Defines the minimum
      allowed orthogonality between vectors in the local space. If a new vector
      does not satisfy this requirement, the space is modified.

    * `"backtrack tolerance`" ``[double]`` **0.** Require a reduction of at
      least this much in the residual norm before accepting a correction.

    * `"backtrack factor`" ``[double]`` **0.5** Multiply the correction by this
      factor each backtracking step.  Note, should be in (0, 1)

    * `"backtrack monitor`" ``[string]`` **monitor either** What norm is
      checked to determine whether backtracking has improved the residual or
      not?  One of `"monitor enorm`", `"monitor L2 residual`", or `'monitor
      either`"

    * `"backtrack max steps`" ``[int]`` **10** Controls how many multiples of
      the backtrack factor are applied before declaring failure.

    * `"backtrack max total steps`" ``[int]`` **1e6** Controls how many total
      backtrack steps may be taken before declaring failure.

    * `"backtrack lag iterations`" ``[int]`` **0** Delay requiring a reduction
      in residual for this many nonlinear iterations.

    * `"backtrack last iteration`" ``[int]`` **1e6** Stop requiring a
      reductiontion in residual after this many nonlinear iterations.

    * `"backtrack fail on bad search direction`" ``[bool]`` **false** If
      backtracking for the full number of "backtrack max steps" is taken, and
      the residual norm has still not be reduced suffiently, this determines
      the behavior.  If true, the solver declares failure.  If false, it takes
      the bad step anyway and hopes to recover in later iterates.

    IF

    * `"Anderson mixing`" ``[bool]`` **false** If true, use Anderson mixing instead of NKA.

    THEN

    * `"relaxation parameter`" ``[double]`` **0.7** The relaxation parameter
      for Anderson mixing.

    END




Solver: NOX
----------------------------------
 Calls Nox nonlinear solvers/JFNK.


Linear Solvers
==============
.. _LinearOperator:
.. _`LinearSolvers`:
.. _`Linear Solver`:
.. _`Inverse`:

 Base class for providing ApplyInverse() methods on operators.

Matrix provides existing Operators with an inverse.  Note this may be
iterative or non-iterative, assembled or non-assembled, approximate or exact to
machine precision.

.. _inverse-typed-spec:
.. admonition:: inverse-typed-spec

   * `"iterative method`" ``[string]`` **optional**
   * `"direct method`" ``[string]`` **optional**
   * `"preconditioning method`" ``[string]`` **optional**

DOCUMENT ME!




Linear Solver: PCG
--------------------
 Preconditioned conjugate gradient method for a linear solver.

.. _iterative-method-pcg-spec:
.. admonition:: iterative-method-pcg-spec

    * `"error tolerance`" ``[double]`` **1.e-6** Tolerance on which to declare success.

    * `"maximum number of iterations`" ``[int]`` **100** Maximum iterations before declaring failure.

    * `"overflow tolerance`" ``[double]`` **3.e50** Error above this value results in failure.

    * `"convergence criterial`" ``[Array(string)]`` **{relative rhs}** A list of
      criteria, any of which can be applied.  Valid include:

      - `"relative rhs`" : measure error relative to the norm of the RHS vector
      - `"relative residual`" : measure error relative to the norm of the residual
      - `"absolute residual`" : measure error directly, norm of error
      - `"make one iteration`" : require at least one iteration to be performed before declaring success




Linear Solver: GMRES
--------------------
 Generalized minimum residual method for a linear solver.

Based on the methods of Yu. Kuznetsov, 1968; Y.Saad, 1986.  Deflated version of
GMRES is due to R.Morgan, GMRES with deflated restarting, 2002 SISC; S.Rollin,
W.Fichtner, Improving accuracy of GMRES with deflated restarting, 2007 SISC.

.. _iterative-method-gmres-spec:
.. admonition:: iterative-method-gmres-spec

    * `"error tolerance`" ``[double]`` **1.e-6** Tolerance on which to declare success.

    * `"maximum number of iterations`" ``[int]`` **100** Maximum iterations before declaring failure.

    * `"overflow tolerance`" ``[double]`` **3.e50** Error above this value results in failure.

    * `"convergence criterial`" ``[Array(string)]`` **{relative rhs}** A list of
      criteria, any of which can be applied.  Valid include:

      - `"relative rhs`" : measure error relative to the norm of the RHS vector
      - `"relative residual`" : measure error relative to the norm of the residual
      - `"absolute residual`" : measure error directly, norm of error
      - `"make one iteration`" : require at least one iteration to be performed before declaring success

    * `"size of Krylov space`" ``[int]`` **10** Size of the Krylov space used to span the residual.

    * `"controller training start`" ``[int]`` **0** Start iteration for determining
      convergence rates. (Add more please!)

    * `"controller training end`" ``[int]`` **3** Start iteration for determining
      convergence rates. (Add more please!)

    * `"preconditioning strategy`" ``[string]`` **left** Valid are "left" and
      "right"-type preconditioning (see Saad 1986)

    * `"maximum size of deflation space`" ``[int]`` **0** Size of the deflation space, see Rollin et al.




Linear Solver: NKA
--------------------
 Uses NKA method as a linear solver.

This is effectively equivalent to GMRES with a rolling restart, where vectors
fall off the end of the space.

.. _iterative-method-nka-spec:
.. admonition:: iterative-method-nka-spec

    * `"error tolerance`" ``[double]`` **1.e-6** Tolerance on which to declare success.

    * `"maximum number of iterations`" ``[int]`` **100** Maximum iterations before declaring failure.

    * `"overflow tolerance`" ``[double]`` **3.e50** Error above this value results in failure.

    * `"convergence criterial`" ``[Array(string)]`` **{relative rhs}** A list of
      criteria, any of which can be applied.  Valid include:

      - `"relative rhs`" : measure error relative to the norm of the RHS vector
      - `"relative residual`" : measure error relative to the norm of the residual
      - `"absolute residual`" : measure error directly, norm of error
      - `"make one iteration`" : require at least one iteration to be performed before declaring success

    * `"max nka vectors`" ``[int]`` **10** Size of the NKA space used to span the residual, conceptually equivalent to the size of the Krylov space.

    * `"nka vector tolerance`" ``[double]`` **0.05** Vectors whose dot product are within this tolerance are considered parallel, and therefore the old vector is thrown out.





Linear Solver: Amesos
---------------------
  Direct solvers via Trilinos.

.. warning:: undocumented




Linear Solver: Amesos
---------------------
  Direct solvers via Trilinos.

.. warning:: undocumented




Linear Solver: Belos (GMRES)
----------------------------
 Trilinos/Belos implementations of iterative methods.
Includes GMRES

.. warning:: undocumented





Preconditioners
===============
.. _Preconditioner:

 A base class for assembled preconditioners.

Identity
--------
 Identity as a preconditioner.

Simply copies the input vector to the output -- uses the Identity matrix as a
preconditioner.

This is provided when using the `"preconditioning method`"=`"identity`" in the
`Preconditioner`_ spec.

No parameters are required.




Diagonal
--------
 Diagonal preconditioner.

Simply applys the pointwise inverse of the diagonal of the matrix as an
extremely cheap matrix.

This is provided when using the `"preconditioning method`"=`"diagonal`" in the
`Preconditioner`_ spec.

No parameters are required.




Block ILU
---------
 Ifpack suite of preconditioners, including block ILU.

The Ifpack (Incomplete Factorization Package) from Trilinos provides
additive-Schwarz-based incomplete factorization methods, including ILU and many
others.

The method is provided in the parent list, e.g. `"preconditioning method`" =
`"ifpack: METHOD`" when calling the Preconditioner factory.

Valid methods include:

- `"point relaxation`" : Additive Schwarz + relaxation
- `"block relaxation`" : Additive Schwarz + block relaxation
- `"Amesos`" : Additive Schwarz with Amesos on block
- `"IC`" : Additive Schwarz + Incomplete Cholesky (symmetry required)
- `"ICT`" : Additive Schwarz + tolerance based Incomplete Cholesky (symmetry required)
- `"ILU`" : Additive Schwarz + Incomplete LU
- `"ILUT`" : Additive Schwarz + tolerance based Incomplete LU
- `"block ilu`" : is the same as `"ILU`", as this is the legacy Amanzi name.

Note that all of these can be used without Additive Schwarz by appending
"stand-alone".

The full list of relevant parameters is somewhat method-dependent, and is documented extensively at

.. url:: https://docs.trilinos.org/dev/packages/ifpack/doc/html/index.html

Here we document a subset of the most frequently used parameters -- advanced
users should read the Ifpack User Guide above to see all options.

.. _preconditioner-ifpack-ilu-spec:
.. admonition:: preconditioner-ifpack-ilu-spec:

    * `"schwarz: combine mode`" ``[string]`` **Add** Note that `"Zero`" may
      perform better for nonsymmetric cases.
    * `"overlap`" ``[int]`` **0** overlap of the Additive Schwarz
    * `"fact: relax value`" ``[double]`` **0.0** If nonzero, dropped values are added to the diagonal (times this factor).
    * `"fact: absolute threshold`" ``[double]`` **0.0** Defines the value to
      add to each diagonal element (times the sign of the actual diagonal
      element).
    * `"fact: relative threshold`" ``[double]`` **1.0** Multiplies the
      diagonal by this value before checking the threshold.
    * `"fact: level-of-fill`" ``[int]`` **0**

.. _preconditioner-ifpack-relaxation-spec:
.. admonition:: preconditioner-ifpack-relaxation-spec:

    * `"schwarz: combine mode`" ``[string]`` **Add** Note that `"Zero`" may
      perform better for nonsymmetric cases.
    * `"overlap`" ``[int]`` **0** overlap of the Additive Schwarz
    * `"relaxation: type`" ``[string]`` **Jacobi**
    * `"relaxation: sweeps`" ``[int]`` **1**
    * `"relaxation: damping factor`" ``[double]`` **1.0**

.. _preconditioner-ifpack-amesos-spec:
.. admonition:: preconditioner-amesos-relaxation-spec:

    * `"schwarz: combine mode`" ``[string]`` **Add** Note that `"Zero`" may
      perform better for nonsymmetric cases.
    * `"overlap`" ``[int]`` **0** overlap of the Additive Schwarz
    * `"amesos: solver type`" ``[string]`` **Amesos_Klu**
    
Example:

.. code-block:: xml

  <ParameterList name="block ilu parameters">
    <Parameter name="fact: relax value" type="double" value="1.0"/>
    <Parameter name="fact: absolute threshold" type="double" value="0.0"/>
    <Parameter name="fact: relative threshold" type="double" value="1.0"/>
    <Parameter name="fact: level-of-fill" type="int" value="0"/>
    <Parameter name="overlap" type="int" value="0"/>
    <Parameter name="schwarz: combine mode" type="string" value="Add"/>
    </ParameterList>
  </ParameterList>




Boomer AMG and Euclid
----------------------
 Hypre based preconditioners include Algebraic MultiGrid and global ILU

Boomer AMG is a HYPRE product consisting of a variety of Algebraic Multigrid
methods.  It is accessed through Ifpack.

This is provided when using the `"preconditioning method`"=`"boomer amg`" or
`"preconditioning method`" = `"hypre: boomer amg`" in the `Preconditioner`_
spec.

.. _preconditioner-boomer-amg-spec:
.. admonition:: preconditioner-boomer-amg-spec:

    * `"tolerance`" ``[double]`` **0.** If is not zero, the preconditioner is dynamic
      and approximate the inverse matrix with the prescribed tolerance (in
      the energy norm ???).

    * `"smoother sweeps`" ``[int]`` **3** defines the number of smoothing loops. Default is 3.

    * `"cycle applications`" ``[int]`` **5** defines the number of V-cycles.

    * `"strong threshold`" ``[double]`` **0.5** defines the number of V-cycles. Default is 5.

    * `"relaxation type`" ``[int]`` **6** defines the smoother to be used. Default is 6
      which specifies a symmetric hybrid Gauss-Seidel / Jacobi hybrid method. TODO: add others!

    * `"coarsen type`" ``[int]`` **0** defines the coarsening strategy to be used. Default is 0
      which specifies a Falgout method. TODO: add others!

    * `"max multigrid levels`" ``[int]`` optionally defined the maximum number of multigrid levels.

    * `"use block indices`" ``[bool]`` **false** If true, uses the `"systems of
      PDEs`" code with blocks given by the SuperMap, or one per DoF per entity
      type.

    * `"number of functions`" ``[int]`` **1** Any value > 1 tells Boomer AMG to
      use the `"systems of PDEs`" code with strided block type.  Note that, to use
      this approach, unknowns must be ordered with DoF fastest varying (i.e. not
      the native Epetra_MultiVector order).  By default, it uses the `"unknown`"
      approach in which each equation is coarsened and interpolated independently.

    * `"nodal strength of connection norm`" ``[int]`` tells AMG to coarsen such
      that each variable has the same coarse grid - sometimes this is more
      "physical" for a particular problem. The value chosen here for nodal
      determines how strength of connection is determined between the coupled
      system.  I suggest setting nodal = 1, which uses a Frobenius norm.  This
      does NOT tell AMG to use nodal relaxation.  Default is 0.

    * `"verbosity`" ``[int]`` **0** prints a summary of run time settings and
      timing information to stdout.  `"1`" prints coarsening info, `"2`" prints
      smoothing info, and `"3`'" prints both.

Example:

.. code-block:: xml

  <ParameterList name="boomer amg parameters">
    <Parameter name="tolerance" type="double" value="0.0"/>
    <Parameter name="smoother sweeps" type="int" value="3"/>
    <Parameter name="cycle applications" type="int" value="5"/>
    <Parameter name="strong threshold" type="double" value="0.5"/>
    <Parameter name="coarsen type" type="int" value="0"/>
    <Parameter name="relaxation type" type="int" value="3"/>
    <Parameter name="verbosity" type="int" value="0"/>
    <Parameter name="number of functions" type="int" value="1"/>
  </ParameterList>


Euclid is a Parallel Incomplete LU, provided as part of the HYPRE project
through the Ifpack interface.

This is provided when using the `"preconditioning method`"=`"euclid`" or
=`"hypre: euclid`" in the `Preconditioner`_ spec.

.. _preconditioner-euclid-spec:
.. admonition:: preconditioner-euclid-spec:

    * `"ilu(k) fill level`" ``[int]`` **1** The factorization level.
    * `"ilut drop tolerance`" ``[double]`` **0** Defines a drop tolerance relative to the largest absolute value of any entry in the row being factored.
    * `"rescale row`" ``[bool]`` **false** If true, values are scaled prior to factorization so that largest value in any row is +1 or -1. Note that this can destroy matrix symmetry.
    * `"verbosity`" ``[int]`` **0** Prints a summary of runtime settings and timing information to stdout.





ML (Trilinos AMG)
-----------------
 Trilinos ML smoothed aggregation multigrid.

This is provided when using the `"preconditioning method`"=`"ml`" in the
`Preconditioner`_ spec.

.. warning:: no input spec defined

See also: https://trilinos.github.io/pdfs/mlguide5.pdf

Example:

.. code-block:: xml

   <ParameterList name="ml parameters">
     <Parameter name="ML output" type="int" value="0"/>
     <Parameter name="aggregation: damping factor" type="double" value="1.33"/>
     <Parameter name="aggregation: nodes per aggregate" type="int" value="3"/>
     <Parameter name="aggregation: threshold" type="double" value="0.0"/>
     <Parameter name="aggregation: type" type="string" value="Uncoupled"/>
     <Parameter name="coarse: type" type="string" value="Amesos-KLU"/>
     <Parameter name="coarse: max size" type="int" value="128"/>
     <Parameter name="coarse: damping factor" type="double" value="1.0"/>
     <Parameter name="cycle applications" type="int" value="2"/>
     <Parameter name="eigen-analysis: iterations" type="int" value="10"/>
     <Parameter name="eigen-analysis: type" type="string" value="cg"/>
     <Parameter name="max levels" type="int" value="40"/>
     <Parameter name="prec type" type="string" value="MGW"/>
     <Parameter name="smoother: damping factor" type="double" value="1.0"/>
     <Parameter name="smoother: pre or post" type="string" value="both"/>
     <Parameter name="smoother: sweeps" type="int" value="2"/>
     <Parameter name="smoother: type" type="string" value="Gauss-Seidel"/>
   </ParameterList>

 



Other Common Specs
##################

IOEvent
=======
 IOEvent: base time/timestep control determing when in time to do something.

The IOEvent is used for multiple objects that need to indicate simulation times or cycles on which to do something.

.. _io-event-spec:
.. admonition:: io-event-spec

   * `"cycles start period stop`" ``[Array(int)]`` **optional**

      The first entry is the start cycle, the second is the cycle
      period, and the third is the stop cycle or -1, in which case there
      is no stop cycle. A visualization dump is written at such
      cycles that satisfy cycle = start + n*period, for n=0,1,2,... and
      cycle < stop if stop != -1.0.

   * `"cycles start period stop 0`" ``[Array(int)]`` **optional** 

      If multiple cycles start period stop parameters are needed, then use these
      parameters.  If one with 0 is found, then one with 1 is looked for, etc,
      until the Nth one is not found.

   * `"cycles`" ``[Array(int)]``  **optional**
  
      An array of discrete cycles that at which a visualization dump is
      written.

   * `"times start period stop`" ``[Array(double)]`` **optional** 

      The first entry is the start time, the second is the time period,
      and the third is the stop time or -1, in which case there is no
      stop time. A visualization dump is written at such times that
      satisfy time = start + n*period, for n=0,1,2,... and time < stop
      if stop != -1.0.

   * `"times start period stop units`" ``[string]`` **s** 

      Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
    
    * `"times start period stop 0`" ``[Array(double)]`` **optional**

      If multiple start period stop parameters are needed, then use this these
      parameters with N=0,1,2.  If one with 0 is found, then one with 1 is
      looked for, etc, until the Nth one is not found.

    * `"times start period stop 0 units`" ``[string]`` **s** 

      Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
      See above for continued integer listings.

    * `"times`" ``[Array(double)]`` **optional** 

      An array of discrete times that at which a visualization dump
      shall be written.

    * `"times units`" ``[string]`` **s** 

      Units corresponding to this spec.  One of `"s`", `"d`", `"yr`", or `"yr 365`"
    
 


Verbose Object
==============
 VerboseObject: a controller for writing log files on multiple cores with varying verbosity.

This allows control of log-file verbosity for a wide variety of objects
and physics.

* `"verbosity level`" ``[string]`` **GLOBAL_VERBOSITY**, `"low`", `"medium`", `"high`", `"extreme`"

   The default is set by the global verbosity spec, (fix me!)  Typically,
   `"low`" prints out minimal information, `"medium`" prints out errors and
   overall high level information, `"high`" prints out basic debugging, and
   `"extreme`" prints out local debugging information.

Note: while there are other options, users should typically not need them.
Instead, developers can use them to control output.
   
Example:

.. code-block:: xml

  <ParameterList name="verbose object">
    <Parameter name="verbosity level" type="string" value="medium"/>
    <Parameter name="name" type="string" value="my header"/>
    <Parameter name="hide line prefix" type="bool" value="false"/>
    <Parameter name="write on rank" type="int" value="0"/>
  </ParameterList>




Debugger
========
 A mesh and vector structure aware utility for printing info.

This is a utility that makes it easier for the user to control output written
to the screen.  It allows the user to provide element IDs, and then provides
functionality for a PK to write mesh geometry information and vector values of
those elements to screen based upon verbosity levels.

Note, most information is only written if the owning object's verbosity level
from the `"Verbose Object`" spec is set to `"high`" or higher.

.. debugger-spec:
.. admonition:: debugger-spec

    * `"debug cells`" ``[Array(int)]`` For each global ID of a cell provided
      here, controls writing of vectors inside of the using PK.

    * `"debug faces`" ``[Array(int)]`` For each global ID of a face provided
      here, writes all adjoining cell information as if each cell was included
      in `"debug cells`".




Residual Debugger
=================


Debugging object for writing vectors to file within an iterative
process for use with vis tools.

.. _residual-debugger-spec:
.. admonition:: residual-debugger-spec

    * `"file name base`" ``[string]`` **amanzi_dbg** Prefix for output filenames.

    INCLUDES:
    - ``[io-event-spec]`` An IOEvent_ spec




   

Function
===================
 Function: base class for all functions of space and time.
Analytic, algabraic functions of space and time are used for a variety of
purposes, including boundary conditions, initial conditions, and independent
variables.

For initial conditions, functions are prescribed of space only, i.e.

:math:`u = f(x,y,z)`

For boundary conditions and independent variables, functions are also a
function of time:

:math:`u = f(t,x,y,z)`

Note, this does not follow the `"typed`" format for legacy reasons.




It is straightforward to add new functions as needed.

Constant Function
-------------------------
 FunctionConstant: Implements the Function interface using a constant value.

Constant function is defined as :math:`f(x) = a`, for all :math:`x`. 

* `"value`" ``[double]`` The constant to be applied.

Example:

.. code-block:: xml

  <ParameterList name="function-constant">
    <Parameter name="value" type="double" value="1.0"/>
  </ParameterList>


  

Tabular Function
-------------------------
 FunctionTabular: Piecewise-defined function.

A piecewise function of one variable.

A tabular function is tabulated on a series of intervals; given values
:math:`{{x_i}}, {{y_i}},, i=0, ... n-1` and functional forms :math:`{{f_j}},,
j=0, ... n-2` a tabular function :math:`f(x)` is defined as:

.. math::
  \begin{matrix}
  f(x) &=& y_0, & x \le x_0,\\
  f(x) &=& f_{{i-1}}(x)  & x \in (x_{{i-1}}, x_i],\\
  f(x) &=& y_{{n-1}}, & x > x_{{n-1}}.
  \end{matrix}

The functional forms :math:`{f_j}` may be constant, which uses the left endpoint, i.e.

:math:`f_i(x) = y_i`,

linear, i.e.

:math:`f_i(x) = ( y_i * (x - x_i) + y_{{i+1}} * (x_{{i+1}} - x) ) / (x_{{i+1}} - x_i)`

or arbitrary, in which the :math:`f_j` must be provided.

The :math:`x_i` and :math:`y_i` may be provided in one of two ways -- explicitly in the input spec or from an HDF5 file.  The length of these must be equal, and the :math:`x_i` must be monotonically increasing.  Forms, as defined on intervals, must be of length equal to the length of the :math:`x_i` less one.

Explicitly specifying the data:

.. _function-tabular-spec:
.. admonition:: function-tabular-spec

   * `"x values`" ``[Array(double)]`` the :math:`x_i`
   * `"y values`" ``[Array(double)]`` the :math:`y_i`

   * `"forms`" ``[Array(string)]`` **optional** Form of the interpolant, either
     `"constant`", `"linear`", or `"USER_DEFINED`" Default is linear for each *
     interval.  Note the length of this array should be one per interval, or
     one less than len of x and y values.

   * `"USER_DEFINED`" ``[function-spec]`` **optional** user-provided functional
     forms on the interval

   * `"x coordinate`" ``[string]`` **t**, `"x`", `"y`", `"z`" defines which
     coordinate direction the :math:`x_i` are formed, defaulting to time.

The below example defines a function that is zero on interval :math:`(-\infty,\,0]`,
linear on interval :math:`(0,\,1]`, constant (`f(x)=1`) on interval :math:`(1,\,2]`,
square root of `t` on interval :math:`(2,\,3]`,
and constant (`f(x)=2`) on interval :math:`(3,\,\infty]`.

Example:

.. code-block:: xml

  <ParameterList name="function-tabular">
    <Parameter name="x values" type="Array(double)" value="{0.0, 1.0, 2.0, 3.0}"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="y values" type="Array(double)" value="{0.0, 1.0, 2.0, 2.0}"/>
    <Parameter name="forms" type="Array(string)" value="{linear, constant, USER_FUNC}"/>

    <ParameterList name="USER_FUNC">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sqrt"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Loading table from file.  (Note that `"USER_DEFINED`" is not an option here,
but could be made so if requested).

.. _function-tabular-fromfile-spec:
.. admonition:: function-tabular-fromfile-spec

   * `"file`" ``[string]`` filename of the HDF5 data
   * `"x header`" ``[string]`` name of the dataset for the :math:`x_i` in the file
   * `"y header`" ``[string]`` name of the dataset for the :math:`y_i` in the file
   * `"forms`" ``[string]`` **optional**, Form of the interpolant, either
     `"constant`" or `"linear`"

The example below would perform linear-interpolation on the intervals provided by data within the hdf5 file `"my_data.h5`".

Example:

.. code-block:: xml

  <ParameterList name="function-tabular">
    <Parameter name="file" type="string" value="my_data.h5"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="x header" type="string" value="/time"/>
    <Parameter name="y header" type="string" value="/data"/>
  </ParameterList>




Smooth step Function
-------------------------
 FunctionSmoothStep: a smoothed discontinuity.

A smooth :math:`C^2` function `f(x)` on interval :math:`[x_0,\,x_1]` is
defined such that `f(x) = y_0` for `x < x0`, `f(x) = y_1` for `x > x_1`, and
monotonically increasing for :math:`x \in [x_0, x_1]` through cubic
interpolation.

* `"x0`" ``[double]`` First fitting point
* `"y0`" ``[double]`` First fitting value
* `"x1`" ``[double]`` Second fitting point
* `"y1`" ``[double]`` Second fitting value

Example:

.. code-block:: xml

  <ParameterList name="function-smooth-step">
    <Parameter name="x0" type="double" value="0.0"/>
    <Parameter name="y0" type="double" value="0.0"/>
    <Parameter name="x1" type="double" value="1.0"/>
    <Parameter name="y1" type="double" value="2.0"/>
  </ParameterList>




Polynomial Function
-------------------------
 FunctionPolynomial: a polynomial

A generic polynomial function is given by the following expression:

.. math::
  f(x) = \sum_{{j=0}}^n c_j (x - x_0)^{{p_j}}

where :math:`c_j` are coefficients of monomials,
:math:`p_j` are integer exponents, and :math:`x_0` is the reference point.

* `"coefficients`" ``[Array(double)]`` c_j polynomial coefficients
* `"exponents`" ``[Array(int)]`` p_j polynomail exponents
* `"reference point`" ``[double]`` x0 to which polynomial argument is normalized.

Example:

.. code-block:: xml

  <ParameterList name="function-polynomial">
    <Parameter name="coefficients" type="Array(double)" value="{1.0, 1.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 4}"/>
    <Parameter name="reference point" type="double" value="0.0"/>
  </ParameterList>


  

Multi-variable linear Function
------------------------------
 FunctionLinear: a multivariate linear function.

A multi-variable linear function is formally defined by
 
.. math::
  f(x) = y_0 + \sum_{{j=0}}^{{n-1}} g_j (x_j - x_{{0,j}}) 

with the constant term "math:`y_0` and gradient :math:`g_0,\, g_1\,..., g_{{n-1}}`.
If the reference point :math:`x_0` is specified, it must have the same
number of values as the gradient.  Otherwise, it defaults to zero.
Note that one of the parameters in a multi-valued linear function can be time.

* `"y0`" ``[double]`` y_0 in f = y0 + g * (x - x0)
* `"gradient`" ``[Array(double)]`` g in f = y0 + g * (x - x0)
* `"x0`" ``[Array(double)]`` x0 in f = y0 + g * (x - x0)

Conditions:

.. code-block:: python

  len(x0) == len(gradient)


Example:

.. code-block:: xml

  <ParameterList name="function-linear">
    <Parameter name="y0" type="double" value="1.0"/>
    <Parameter name="gradient" type="Array(double)" value="{{1.0, 2.0, 3.0}}"/>
    <Parameter name="x0" type="Array(double)" value="{{2.0, 3.0, 1.0}}"/>
  </ParameterList>


  

Separable Function
------------------
 FunctionSeparable: f(x,y) = f1(x)*f2(y)

A separable function is defined as the product of other functions such as

.. math::
  f(x_0, x_1,...,x_{{n-1}}) = f_1(x_0)\, f_2(x_1,...,x_{{n-1}})

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist.

* `"function1`" ``[function-spec]`` f_1 in f(x) = f_1(x0) * f_2(x1...)
* `"function2`" ``[function-spec]`` f_2 in f(x) = f_1(x0) * f_2(x1...)


.. code-block:: xml

  <ParameterList name="function-separable">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>




Additive Function
------------------
 FunctionAdditive: f(x,y) = f1(x,y) + f2(x,y)

An additive function simply adds two other function results together.

.. math::
  f(x) = f_1(x) + f_2(x)

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist.

* `"function1`" ``[function-spec]`` f_1 in f(x) = f_1(x) + f_2(x)
* `"function2`" ``[function-spec]`` f_2 in f(x) = f_1(x) + f_2(x)

Example:

.. code-block:: xml

  <ParameterList name="function-additive">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>



Multiplicative Function
--------------------------
 FunctionMultiplicative: f(x,y) = f1(x,y) * f2(x,y)

A multiplicative function simply multiplies two other function results together.

.. math::
  f(x) = f_1(x) * f_2(x)

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist.

* `"function1`" ``[function-spec]`` f_1 in f(x) = f_1(x) + f_2(x)
* `"function2`" ``[function-spec]`` f_2 in f(x) = f_1(x) + f_2(x)

Example:

.. code-block:: xml

  <ParameterList name="function-multiplicative">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>



Composition Function
--------------------------
 FunctionComposition: f(x,y) = f1(x,y) * f2(x,y)

Function composition simply applies one function to the result of another.

.. math::
  f(x) = f_1( f_2(x) )

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist.

* `"function1`" ``[function-spec]`` f_1 in f(x) = f_1(f_2(x))
* `"function2`" ``[function-spec]`` f_2 in f(x) = f_1(f_2(x))


.. code-block:: xml

  <ParameterList name="function-composition">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>



Piecewise Bilinear Function
---------------------------
 FunctionBilinear: a piecewise bilinear function.

A piecewise bilinear function extends the linear form of the tabular function to two variables.

Define :math:`i(x) = i : x_i < x <= x_{{i+1}}` and similarly :math:`j(y) = j : y_j < y <= y_{{j+1}}` for monotonically increasing :math:`x_i` and :math:`y_j`.

Given a two-dimensional array :math:`u_{{i,j}}`, :math:`f` is then defined by
bilinear interpolation on :math:`u_{{i(x),j(y)}}, u_{{i(x)+1,j(y)}},
u_{{i(x),j(y)+1}}, u_{{i(x)+1,j(y)+1}}, if :math:`(x,y)` is in
:math:`[x_0,x_n] \times [y_0,y_m]`, linear interpolation if one of :math:`x,y`
are out of those bounds, and constant at the corner value if both are out of
bounds.
 
* `"file`" ``[string]`` HDF5 filename of the data
* `"row header`" ``[string]`` name of the row dataset, the :math:`x_i`
* `"row coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
* `"column header`" ``[string]`` name of the column dataset, the :math:`y_i`
* `"column coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
* `"value header`" ``[string]`` name of the values dataset, the :math:`u_{{i,j}}`

Example:

.. code-block:: xml

  <ParameterList name="function-bilinear">
    <Parameter name="file" type="string" value="pressure.h5"/>
    <Parameter name="row header" type="string" value="/time"/>
    <Parameter name="row coordinate" type="string" value="t"/>
    <Parameter name="column header" type="string" value="/x"/>
    <Parameter name="column coordinate" type="string" value="x"/>
    <Parameter name="value header" type="string" value="/pressure"/>
  </ParameterList>




Distance Function
-----------------
 FunctionDistance: distance from a reference point.

A distance function calculates distance from reference point :math:`x_0`
using by the following expression:

.. math::
  f(x) = \sqrt( \sum_{j=0}^{n} m_j (x_j - x_{0,j})^2 )

Note that the first parameter in :math:`x` can be time.

* `"x0`" ``[Array(double)]`` Point from which distance is measured.
* `"metric`" ``[Array(double)]`` Linear scaling metric, typically all 1s.

Here is an example of a distance function using isotropic metric:

Example:
.. code-block:: xml

  <ParameterList name="function-distance">
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="metric" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
  </ParameterList>




Monomial Function
-----------------
 FunctionMonomial: a multivariate monomial function.

A multi-variable monomial function is given by the following expression:

.. math::
  f(x) = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}

with the constant factor :math:`c`, the reference point :math:`x_0`, and
integer exponents :math:`p_j`. 
Note that the first parameter in :math:`x` can be time.

* `"c`" ``[double]`` c in f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}
* `"x0`" ``[Array(double)]`` x0 in f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}
* `"exponents`" ``[Array(int)]`` p in f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}

Conditions:

.. code-block:: python

  len(x0) == len(exponents)

Here is an example of monomial of degree 6 in three variables:

.. code-block:: xml

  <ParameterList name="function-monomial">
    <Parameter name="c" type="double" value="1.0"/>
    <Parameter name="x0" type="Array(double)" value="{1.0, 3.0, 0.0}"/>
    <Parameter name="exponents" type="Array(int)" value="{2, 3, 1}"/>
  </ParameterList>




Standard Math Function
----------------------
 FunctionStandardMath: provides access to many common mathematical functions.
These functions allow to set up non-trivial time-dependent boundary conditions 
which increases a set of analytic solutions that can be used in convergence 
analysis tests.

.. math::
  f(x) = A * operator( p * (x - s) )

or

.. math::
  f(x) = A * operator(x-s, p)

Note that these operate only on the first coordinate, which is often time.
Function composition can be used to apply these to other coordinates (or
better yet a dimension could/should be added upon request).

* `"operator`" ``[string]`` specifies the name of a standard mathematical function.
  Available options are `"cos`", `"sin`", `"tan`", `"acos`", `"asin`", `"atan`", 
  `"cosh`", `"sinh`", `"tanh`", `"exp`", `"log`", `"log10`", `"sqrt`", `"ceil`",
  `"fabs`", `"floor`", `"mod`", and `"pow`".

* `"amplitude`" ``[double]`` specifies a multiplication factor `a` in formula `a f(x)`. 
  The multiplication factor is ignored by function `mod`. Default value is 1.

* `"parameter`" ``[double]`` **1.0** specifies additional parameter `p` for
  math functions with two arguments. These functions are `"a pow(x[0], p)`"
  and `"a mod(x[0], p)`".  Alternative, scales the argument before
  application, for use in changing the period of trig functions.

* `"shift`" ``[double]`` specifies a shift of the function argument. Default is 0.

Example:

.. code-block:: xml

  <ParameterList name="function-standard-math">
    <Parameter name="operator" type="string" value="sqrt"/>
    <Parameter name="amplitude" type="double" value="1e-7"/>
    <Parameter name="shift" type="double" value="0.1"/>
  </ParameterList>

This example defines function `1e-7 sqrt(t-0.1)`.
 




Operator
========

 Operator represents a linear map, and typically encapsulates a discretization.
``Operator`` represents a map from linear space X to linear space Y.  Typically,
this map is a linear map, and encapsulates much of the discretization involved
in moving from continuous to discrete equations. The spaces X and Y are described
by CompositeVectors (CV). A few maps X->Y are supported.

An ``Operator`` provides an interface for applying both the forward and inverse
linear map (assuming the map is invertible).

Typically the ``Operator`` is never seen by the user; instead the user provides
input information for helper classes based on the continuous mathematical
operator and the desired discretization.  These helpers build the needed
``Operator``, which may include information from multiple helpers (i.e. in the
case of Jacobian Operators for a PDE).

However, one option may be provided by the user, which is related to dealing
with nearly singular operators:

* `"diagonal shift`" ``[double]`` **0.0** Adds a scalar shift to the diagonal
  of the ``Operator``, which can be useful if the ``Operator`` is singular or
  near-singular.



PDE_Accumulation
----------------

``PDE_Accumulation`` assembles the discrete form of :math:`\frac{\partial A}{\partial t}`.

This class is usually used as part of a preconditioner, providing the linearization:

.. math::
  \frac{\partial}{\partial A} \left[ \frac{\partial A}{\partial t} \right]_{A_0} i
  = \frac{|\Omega_E|}{\Delta t}

for a grid element :math:`\Omega_E`.


.. _pde-accumulation-spec:
.. admonition:: pde-accumulation-spec

  * `"entity kind`" ``[string]`` **optional** Typically set by the PK
  * `"number of vectors`" ``[int]`` **optional** Typically set by the PK




PDE_Diffusion
-------------


``PDE_Diffusion`` forms local ``Op`` s and global ``Operator`` s for elliptic equations:

.. math::
  \nabla \cdot k \nabla u

with a variety of discretizations. Note also, for reasons that are one part historical 
and potentially not that valid, this also supports and implementation with an advective 
source, i.e.:

.. math::
  \nabla \cdot K k (\nabla u + b g z)

for gravitational terms in Richards equations.

The input spec for a diffusion operator consists of:

* `"discretization primary`" ``[string]`` See below for supported options.

  - `"fv: default`" the standard two-point flux finite volume discretization
  - `"nlfv: default`" the nonlinear finite volume method of ???
  - MFD methods, including:

    - `"mfd: default`"
    - `"mfd: monotone for hex`"
    - `"mfd: optimized for monotonicity`"
    - `"mfd: two-point flux approximation`"
    - `"mfd: optimized for sparsity`"
    - `"mfd: support operator`"

 Note that the most commonly used are `"fv: default`" for simple test
 problems (this method is not particularly accurate for distorted
 meshes), `"mfd: optimized for sparsity`" for most real problems on
 unstructured meshes, and `"mfd: optimized for monotonicity`" for
 orthogonal meshes with diagonal tensor/scalar coefficients.

* `"gravity`" ``[bool]`` **false** specifies if the gravitational flow term is included

* `"Newton correction`" ``[string]`` specifies a model for non-physical terms 
  that must be added to the matrix. These terms represent Jacobian and are needed 
  for the preconditioner. Available options are `"true Jacobian`" and `"approximate Jacobian`".
  The FV scheme accepts only the first options. The other schemes accept only the second option.

* `"scaled constraint equation`" ``[bool]`` **false** rescales flux continuity equations
  on mesh faces.  These equations are formed without the nonlinear
  coefficient. This option allows us to treat the case of zero nonlinear
  coefficient, which otherwise generates zero rows in the operator, which is
  then singular.  At moment this feature does not work with non-zero gravity
  term.

* `"constraint equation scaling cutoff`" ``[double]`` specifies the cutoff value for
  applying rescaling strategy described above.



 Diffusion generates local Ops and global Operators for an elliptic operator.

A diffusion PDE for generating both local and global operators.

* `"discretization primary`" ``[string]``



Example:

.. code-block:: xml

    <ParameterList name="OPERATOR_NAME">
      <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
      <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
      <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
      <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
      <Parameter name="gravity" type="bool" value="true"/>
      <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
      <Parameter name="nonlinear coefficient" type="string" value="upwind: face"/>
      <Parameter name="Newton correction" type="string" value="true Jacobian"/>
    </ParameterList>




Additional options available only for the MFD family of discretizations include:
  
* `"nonlinear coefficient`" ``[string]`` specifies a method for treating nonlinear
  diffusion coefficient, if any. Available options are `"none`", `"upwind:
  face`", `"divk: cell-face`" (default), `"divk: face`", `"standard: cell`",
  `"divk: cell-face-twin`" and `"divk: cell-grad-face-twin`".  Symmetry
  preserving methods are the divk-family of methods and the classical
  cell-centered method (`"standard: cell`"). The first part of the name
  indicates the base scheme.  The second part (after the semi-column)
  indicates required components of the composite vector that must be provided
  by a physical PK.

* `"discretization secondary`" ``[string]`` specifies the most robust
  discretization method that is used when the primary selection fails to
  satisfy all a priori conditions.  This is typically `"mfd: default`", and is
  used only when an MFD `"discretization primary`" is used.

* `"schema`" ``[Array(string)]`` defines the operator stencil. It is a collection of 
  geometric objects.  Typically this is set by the implementation and is not provided.

* `"preconditioner schema`" ``[Array(string)]`` **{face,cell}** Defines the
  preconditioner stencil.  It is needed only when the default assembling
  procedure is not desirable. If skipped, the `"schema`" is used instead.
  In addition to the default, **{face}** may be used, which forms the Schur
  complement.
   
* `"consistent faces`" ``[list]`` may contain a `"preconditioner`" and
  `"linear operator`" list (see sections Preconditioners_ and LinearSolvers_
  respectively).  If these lists are provided, and the `"discretization
  primary`" is of type `"mfd: *`", then the diffusion method
  UpdateConsistentFaces() can be used.  This method, given a set of cell
  values, determines the faces constraints that satisfy the constraint
  equation in MFD by assembling and inverting the face-only system.  This is
  not currently used by any Amanzi PKs.

* `"diffusion tensor`" ``[string]`` allows us to solve problems with symmetric and 
  non-symmetric (but positive definite) tensors. Available options are *symmetric* 
  (default) and *nonsymmetric*.




Additional options for MFD with the gravity term include:
  
* `"gravity term discretization`" ``[string]`` selects a model for discretizing the 
  gravity term. Available options are `"hydraulic head`" (default) and `"finite volume`". 
  The first option starts with equation for the shifted solution, i.e. the hydraulic 
  head, and derives gravity discretization by the reserve shifting.
  The second option is based on the divergence formula.








PDE_Advection
-------------




``PDE_AdvectionUpwind`` assembles the discrete form of:

.. math::
  \nabla \cdot (q C)

which advects quantity :math:`C` with fluxes :math:`q`.

This is a simple, first-order donor-upwind scheme, and is recommended
for use in diffusion-dominated advection-diffusion equations.





