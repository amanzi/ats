ATS Native XML Input Specification V-dev
****************************************


.. contents:: **Table of Contents**
  :local:
  :depth: 2

  
About the Specification
#######################

ATS, and Amanzi's "native" specificiation, is an xml file following
Trilinos's Teuchos ParameterList schema.  There are only two types of
tags used -- `"Parameter`" and `"ParameterList`".  `"Parameter`"
elements consist of `"name`", `"type`", and `"value`" attributes.
`"ParameterList`" elements use the `"name`" attribute and include
subelements that are other `"ParameterList`" and `"Parameter`"
elements.

The top-most, `"main`" list is read by the code and used to provide
all information needed to run the simulation.  This input spec is
designed for the code, not necessarily for the user.  In general,
avoid writing input files from scratch, and prefer to modify existing
demos or examples.

Here we document the input spec by defining what each possible element
used by the code needs to be well posed.

Specs
=====

In many cases, the input specifies data for a particular parameterized
model, and ATS supports a number of parameterizations.  For example,
initial data might be uniform (the value is required), or linear in y
(the value and its gradient are required).  Where ATS supports a
number of parameterized models for quantity Z, the available models
will be listed by name, and then will be described in the subsequent
section.  For example, the specification for an `"X`" list might look
like:

.. _X-spec:
.. admonition:: X-spec

  * `"Y`" ``[string]`` **default_value** Documentation desribing Y.
  * `"Z`" ``[Z-spec]`` Model for Z, One of `"z1`" or `"z2`" (see below) 

Here, an `"X`" is defined by a `"Y`" and a `"Z`".  The `"Y`" is a
string parameter but the `"Z`" is given by a model (which will require
its own set of parameters).  The options for `"Z`" will then be
described seperately as a `"Z-spec`"


An example of using such a specification:

.. code-block:: xml

    <ParameterList name="X">
      <Parameter name="Y" type="string" value="hello"/>
      <ParameterList name="z2">
        <Parameter name="z2a" type="double" value="0.7"/>
        <Parameter name="z2b" type="int" value="3"/>
      </ParameterList>   
    </ParameterList>   
 

Syntax
======

* Reserved keywords and labels are `"quoted and italicized`" -- these
  labels or values of parameters in user-generated input files must
  match (using XML matching rules) the specified or allowable values.

* User-defined labels are indicated with ALL-CAPS, and are meant to
  represent a typical or default name given by a user - these can be
  names or numbers or whatever serves best the organization of the
  user input data.  Things liked PRESSURE or SURFACE-PONDED_DEPTH can
  be renamed from their defaults if it makes sense to the problem.

* Bold values are default values, and are used if the Parameter
  is not provided.

Naming
======

Variables are named according to a very strong convention.  While some
variables may be overridden by the user, users should choose to follow
these conventions or things like visualization scripts may not behave
as expected.

A variable name looks like one of:

- SUFFIX
- DOMAIN-SUFFIX
- DOMAIN_SET:ID-VARNAME

where:

- When DOMAIN is supplied, it is the "default" mesh, called `"domain`"
  in the mesh list, and otherwise is the name of the mesh (e.g. `"surface`").
- DOMAIN_SET:ID is itself a DOMAIN, where the set defines the
  collection as a whole (from the mesh list) and the ID is defined by
  an index across the collection, e.g. `"column:4`"

Tags indicate the use of a variable at a specific time in the
discretized time interval.  Default tags include `"current`" and
`"next`" indicating the variable at the beginning and end of the
interval, respectively.  Often subcycling and other schemes will
designate special-purpose tags which may be used internally by a
subset of the equations begin solved.  Tags are combined with
variables to indicate a specific data structure,
e.g. `"surface-pressure@NEXT`".

Lastly, derivatives are named using the `"d`" and the `"|`" character,
e.g. `"dsurface-water_content|dsurface-pressure`" is the derivative of
the `"water_content`" variable on the `"surface`" domain with respect
to the `"pressure`" on the same domain.

As a result of these conventions, none of the above individual strings,
(suffixes, domains, domain sets, or IDs) can contain any of the
following reserved characters: `:`, `-`, `|`, `@`.
  

Symbol Index
============

.. include:: symbol_table.rst
  
Main
####


ATS's top level driver is provided the entire input spec as a single list
called `"main`".  That list contains the following required elements:

.. _main-spec:
.. admonition:: main-spec

    * `"cycle driver`" ``[coordinator-spec]``  See below.
    * `"mesh`" ``[mesh-typed-spec-list]`` A list of Mesh_ objects.
    * `"regions`" ``[region-typedinline-spec-list]`` A list of Region_ objects.
    * `"visualization`" ``[visualization-spec-list]`` A list of Visualization_ objects.
    * `"observations`" ``[observation-spec-list]`` An list of Observation_ objects.
    * `"checkpoint`" ``[checkpoint-spec]`` A Checkpoint_ spec.
    * `"PKs`" ``[pk-typedinline-spec-list]`` A list of `Process Kernels`_.
    * `"state`" ``[state-spec]`` A State_ spec.


Coordinator
############

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
         <Parameter name="number of cells" type="Array(int)" value="{100, 1, 100}"/>
         <Parameter name="domain low coordinate" type="Array(double)" value="{0.0, 0.0, 0.0}" />
         <Parameter name="domain high coordinate" type="Array(double)" value="{100.0, 1.0, 10.0}" />
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
   * `"entity LID`" ``[int]`` Local ID of the surface cell that is the top of the column.
   * `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
   * `"deformable mesh`" ``[bool]`` **false**  Used for deformation PKs to allow non-const access.

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
     empty lists.

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

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

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
    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />
    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>
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
    * `"delimiter`" ``[string]`` **COMMA** Delimiter to split columns of the file
    *  `"write interval`" ``[int]`` **1** Interval of observations on which to flush IO files.
    * `"time units`" ``[string]`` **s** Controls the unit of the time column in the observation file.
    * `"domain`" ``[string]`` **"domain"** Can be used to set the communicator which writes (defaults to the standard subsurface domain).
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
involves a field quantity, a reduction reduction operator, a region from which
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

    * `"reduction`" ``[string]`` the type of function to apply to the variable
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

    * `"modifier`" ``[function-typedinline-spec]`` **optional** If provided, defines a
      function used to modify the `"variable`" prior to applying the
      `"reduction`".
        
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





Process Kernels
###############
.. _`PK`:

Process Kernels, or PKs, are the fundamental unit of a model, and
represent a single or system of Partial Differential Equations (PDEs)
or Differential Algebraic Equations (DAEs).  PKs are broadly split
into individual equations (Physical PKs) and systems of equations,
called Multi-Process Coordinators (MPCs).

The PK tree forms the fundamental definition of the entire system of
equations to be solved by the simulator, and is represented by a
single PK or a single MPC which couples other MPCs and/or Physical
PKs.

.. contents:: **List of PKs**
   :local:

Base PKs
========
There are several types of PKs, and each PK has its own valid input
spec.  However, there are three main types of PKs, from which nearly
all PKs derive.  Note that none of these are true PKs and cannot stand
alone.

PK base class
-------------
 The interface for a Process Kernel, an equation or system of equations.

A process kernel represents a single or system of partial/ordinary
differential equation(s) or conservation law(s), and is used as the
fundamental unit for coupling strategies.

Implementations of this interface typically are either an MPC
(multi-process coupler) whose job is to heirarchically couple several
other PKs and represent the system of equations, or a Physical PK,
which represents a single equation.

Note there are two PK specs -- the first is the "typed" spec, which appears in
the "cycle driver" list in the PK tree and has no other parameters other than
its type and its children.  The second is the spec for the base class PK, which
is inherited and included by each actual PK, lives in the "PKs" sublist of
"main", and has all needed parameters.

.. _pk-typed-spec:
.. admonition:: pk-typed-spec

    * `"PK type`" ``[string]`` One of the registered PK types

Example:

.. code-block:: xml

  <ParameterList name="PK tree">
    <ParameterList name="Top level MPC">
      <Parameter name="PK type" type="string" value="strong MPC"/>
      <ParameterList name="sub PK 1">
        ...
      </ParameterList>
      <ParameterList name="sub PK 2">
        ...
      </ParameterList>
      ...
    </ParameterList>
  </ParameterList>


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

    * `"initial time step [s]`" ``[double]`` **1.** Initial time step size [s]

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

   * `"boundary conditions`" ``[list]`` Defaults to Neuman, 0 normal
     flux.  See `Flow-specific Boundary Conditions`_

   * `"permeability type`" ``[string]`` **scalar** This controls the
     number of values needed to specify the absolute permeability.
     One of:

     - `"scalar`" Requires one scalar value.
     - `"horizontal and vertical`" Requires two values, horizontal
       then vertical.
     - `"diagonal tensor`" Requires dim values: {xx, yy} or {xx, yy,
       zz}
     - `"full tensor`". (Note symmetry is required.)  Either {xx, yy,
       xy} or {xx,yy,zz,xy,xz,yz}.

   * `"water retention evaluator`" ``[wrm-evaluator-spec]`` The water
     retention curve.  This needs to go away, and should get moved to
     State.

   IF

   * `"source term`" ``[bool]`` **false** Is there a source term?

   THEN

   * `"source key`" ``[string]`` **DOMAIN-water_source** Typically not
     set, as the default is good. ``[mol s^-1]``
   * `"source term is differentiable`" ``[bool]`` **true** Can the
     source term be differentiated with respect to the primary
     variable?
   * `"explicit source term`" ``[bool]`` **false** Apply the source
     term from the previous time step.

   END

   Math and solver algorithm options:

   * `"diffusion`" ``[pde-diffusion-spec]`` The (forward) diffusion
     operator, see PDE_Diffusion_.

   * `"diffusion preconditioner`" ``[pde-diffusion-spec]``
     **optional** The inverse of the diffusion operator.  See
     PDE_Diffusion_.  Typically this is only needed to set Jacobian
     options, as all others probably should match those in
     `"diffusion`", and default to those values.

   * `"surface rel perm strategy`" ``[string]`` **none** Approach for
     specifying the relative permeabiilty on the surface face.
     `"clobber`" is frequently used for cases where a surface rel
     perm will be provided.  One of:

     - `"none`" : use the upwind direction to determine whether to
       use the boundary face or internal cell
     - `"clobber`" : always use the boundary face rel perm
     - `"max`" : use the max of the boundary face and internal cell
       values
     - `"unsaturated`" : Uses the boundary face when the internal
       cell is not saturated.

   * `"relative permeability method`" ``[string]`` **upwind with Darcy
     flux** Relative permeability is defined on cells, but must be
     calculated on faces to multiply a flux.  There are several
     methods commonly used.  Note these can significantly change
     answers -- you don't want to change these unless you know what
     they mean.  One of:

     - `"upwind with Darcy flux`" First-order upwind method that is
       most common
     - `"upwind with gravity`" Upwinds according to the gravitational
       flux direction
     - `"cell centered`" This corresponds to the harmonic mean, and is
       most accurate if the problem is always wet, but has issues
       when it is dry.
     - `"arithmetic mean`" Face value is the mean of the neighboring
       cells.  Not a good method.

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

   * `"limit correction to pressure change when crossing atmospheric [Pa]`" ``[double]`` **-1**
     If > 0, this limits an iterate's max pressure change
     to this value when they cross atmospheric pressure.  Not usually helpful.

   Discretization / operators / solver controls:

   * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` **optional**
     The inverse of the accumulation operator.  See PDE_Accumulation_.
     Typically not provided by users, as defaults are correct.

   * `"absolute error tolerance`" ``[double]`` **2750.0** in units of [mol].

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

    * `"limit correction to pressure change when crossing atmospheric [Pa]`" ``[double]`` **-1**
      If > 0, this limits an iterate's max pressure change
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



The advection-diffusion equation for component *i* in partially saturated porous media may be written as

.. math::
  \frac{\partial (\phi s_l C_i)}{\partial t}
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q} C_i)
  + \boldsymbol{\nabla} \cdot (\phi s_l\, (\boldsymbol{D^*}_l + \tau \boldsymbol{D}_i) \boldsymbol{\nabla} C_i) + Q_s,

The advection-diffusion equation for component *i* in the surface may be written as

.. math::
  \frac{\partial (C_i)}{\partial t}
  =
  - \boldsymbol{\nabla} \cdot (\boldsymbol{q_s} C_i)
  + \boldsymbol{\nabla} \cdot ( (\boldsymbol{D^*}_l + \tau \boldsymbol{D}_i) \boldsymbol{\nabla} C_i) + Q_s,

.. _transport-spec:
.. admonition:: transport-spec

    * `"PK type`" ``[string]`` **"transport ats"**

    * `"domain name`" ``[string]`` **domain** specifies mesh name that defines
      the domain of this PK.

    * `"component names`" ``[Array(string)]`` No default. Provides the names of the
      components that will be transported. Must be in the order: aqueous, gaseous, solid.

    * `"number of aqueous components`" ``[int]`` **-1** The total number of
      aqueous components.  Default value is the length of `"component names`"

    * `"number of gaseous components`" ``[int]`` **0** The total number of
      gaseous components.

    * `"boundary conditions`" ``[transport-bc-spec]`` Boundary conditions for
      transport are dependent on the boundary conditions for flow. See
      `Flow-specific Boundary Conditions`_ and `Transport-specific Boundary Conditions`_

    * `"component molar masses`" ``[Array(double)]`` No default. Molar mass of
      each component.

    * `"molecular diffusion`" ``[molecular-diffusion-spec]`` defines names of
      solutes in aqueous and gaseous phases and related diffusivity values.

    * "material properties" [material-properties-spec-list] Defines material
      properties see below).

    Source terms:

    * `"source terms`" [transport-source-spec-list] Provides solute source.

    Physical model and assumptions:

    * `"physical models and assumptions`" [material-properties-spec] Defines material properties.

    * `"effective transport porosity`" [bool] If *true*, effective transport porosity
      will be used by dispersive-diffusive fluxes instead of total porosity.
      Default is *false*.

    Math and solver algorithm options:

    * `"diffusion`" ``[pde-diffusion-spec]`` Diffusion drives the distribution.
      Typically we use finite volume here.  See PDE_Diffusion_

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` Inverse of the
      above.  Likely only Jacobian term options are needed here, as the others
      default to the same as the `"diffusion`" list.  See PDE_Diffusion_.

    * `"inverse`" ``[inverse-typed-spec]`` Inverse_ method for the solve.

    * `"cfl`" [double] Time step limiter, a number less than 1. Default value is 1.

    * `"spatial discretization order`" [int] defines accuracy of spatial discretization.
      It permits values 1 or 2. Default value is 1.

    * `"temporal discretization order`" [int] defines accuracy of temporal discretization.
      It permits values 1 or 2 and values 3 or 4 when expert parameter
      `"generic RK implementation`" is set to true. Note that RK3 is not monotone.
      Default value is 1.

    * `"reconstruction`" [list] collects reconstruction parameters. The available options are
      describe in the separate section below.

    * `"transport subcycling`" ``[bool]`` **true** The code will default to
      subcycling for transport within the master PK if there is one.


    Developer parameters:

    * `"enable internal tests`" [bool] turns on various internal tests during
      run time. Default value is `"false`".

    * `"generic RK implementation`" [bool] leads to generic implementation of
      all Runge-Kutta methods. Default value is `"false`".

    * `"internal tests tolerance`" [double] tolerance for internal tests such as the
      divergence-free condition. The default value is 1e-6.

    * `"runtime diagnostics: solute names`" [Array(string)] defines solutes that will be
      tracked closely each time step if verbosity `"high`". Default value is the first
      solute in the global list of `"aqueous names`" and the first gas in the global list
      of `"gaseous names`".

    * `"runtime diagnostics: regions`" [Array(string)] defines a boundary region for
      tracking solutes. Default value is a seepage face boundary, see Flow PK.

    KEYS

    - `"saturation liquid`" This variable is a multiplier in in the
      accumulation term. For subsurface transport, this will typically be the
      saturation (`"saturation_liquid`"). For surface transport, this will
      typically be the ponded depth (`"ponded_depth`").

    - `"previous saturation liquid`"

    - `"molar density liquid`"  Transport is solved
      for concentrations in units of mol fractions. Molar density is needed for conversion.

    - `"water flux`"

    - `"water source`" Defines the water injection rate [mol H2O m^-2 s^-1] in
      surface and [mol H2O m^-3 s^-1] in subsurface) which applies to
      concentrations specified by the `"geochemical conditions`".  Note that if
      this PK is coupled to a surface flow PK, the unit of the water source
      there *must* be in [mol m^-2 s^-1], *not* in [m s^-1] as is an option for
      that PK (e.g. `"water source in meters`" must be set to `"false`" in the
      overland flow PK).

      The injection rate of a solute [molC s^-1], when given as the product of
      a concentration and a water source, is evaluated as:

      Concentration [mol C L^-1] *
        1000 [L m^-3] of water *
        water source [mol H2O m^-3 s^-1] *
        volume of injection domain [m^3] /
        molar density of water [mol H2O m^-3]


.. _molecular-diffusion-spec:
.. admonition:: molecular-diffusion-spec

   * `"aqueous names`" ``[Array(string)]`` List of aqueous component names to
     be diffused.
   * `"aqueous values`" ``[Array(string)]`` Diffusivities of each component.


.. code-block:: xml

   <ParameterList name="molecular diffusion">
     <Parameter name="aqueous names" type=Array(string)" value="{CO2(l),Tc-99}"/>
     <Parameter name="aqueous values" type=Array(double)" value="{1e-8,1e-9}"/>
   </ParameterList>


.. _material-properties-spec:
.. admonition:: material-properties-spec

   * `"region`" ``[Array(string)]`` Defines geometric regions for material SOIL.

   * `"model`" ``[string]`` **scalar** Defines dispersivity model.  One of:

     - `"scalar`" : scalar dispersivity
     - `"Bear`" : dispersion split into along- and across- velocity
     - `"Burnett-Frind`"
     - `"Lichtner-Kelkar-Robinson`"

   * `"parameters for MODEL`" ``[list]`` where `"MODEL`" is the model name.

   IF model == scalar

   ONE OF

   * `"alpha`" ``[double]`` defines dispersivity in all directions, [m].

   OR

   * `"dispersion coefficient`" ``[double]`` defines dispersion coefficient [m^2 s^-1].

   END

   ELSE IF model == Bear

   * `"alpha_l`" ``[double]`` defines dispersion in the direction of Darcy velocity, [m].
   * `"alpha_t`" ``[double]`` defines dispersion in the orthogonal direction, [m].

   ELSE IF model == Burnett-Frind

   * `"alphaL`" ``[double]`` defines the longitudinal dispersion in the direction
     of Darcy velocity, [m].
   * `"alpha_th`" ``[double]`` Defines the transverse dispersion in the horizonal
     direction orthogonal directions, [m].
   * `"alpha_tv`" ``[double]`` Defines dispersion in the orthogonal directions,
     [m].  When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in
     the direction of the Darcy velocity.

   ELSE IF model == Lichtner-Kelker-Robinson

   * `"alpha_lh`" ``[double]`` defines the longitudinal dispersion in the
     horizontal direction, [m].
   * `"alpha_lv`" ``[double]`` Defines the longitudinal dispersion in the vertical
     direction, [m].  When `"alpha_lh`" equals to `"alpha_lv`", we obtain
     dispersion in the direction of the Darcy velocity.
   * `"alpha_th`" ``[double]`` Defines the transverse dispersion in the horizontal
     direction orthogonal directions, [m].
   * `"alpha_tv" ``[double]`` Defines dispersion in the orthogonal directions.
     When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in the
     direction of the Darcy velocity.

   END

   * `"aqueous tortuosity`" ``[double]`` Defines tortuosity for calculating
     diffusivity of liquid solutes, [-].

   * `"gaseous tortuosity`" ``[double]`` Defines tortuosity for calculating
     diffusivity of gas solutes, [-].


.. _transport-source-spec:
.. admonition:: transport-source-spec

   * `"component mass source`" ``[list]``  Defines solute source injection rate.

     * `"spatial distribution method`" ``[string]`` One of:

        - `"volume`", source is considered as extensive quantity [molC s^-1] and is evenly distributed across the region.
        - `"none`", source is considered as intensive quantity. [molC m^-2 s^-1] in surface and [molC m^-3 s^-1] in subsurface

     * `"geochemical`" ``[list]``  Defines a source by setting solute concentration for all components (in moles/L) and an injection
       rate given by the water source.  Currently, this option is only available for Alquimia provided geochemical conditions.

       - `"geochemical conditions`" ``[Array(string)]`` List of geochemical constraints providing concentration for solute injection.




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
    - `"flux`" **DOMAIN-water_flux** The water flux :math:`\mathbf{q}` used in advection. `[mol s^-1]`
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

    * `"conserved quantity key`" ``[string]`` **DOMAIN-snow_water_equivalent**
      Sets the default conserved quantity key, so this is likely not supplied
      by the user. `[m]`
    * `"snow density key`" ``[string]`` **DOMAIN-density** Default snow density
      key. `[kg m^-3]`
    * `"snow age key`" ``[string]`` **DOMAIN-age** Default snow age key. `[d]`
    * `"new snow key`" ``[string]`` **DOMAIN-source** Default new snow key. `[m SWE s^-1]`
    * `"area fractions key`" ``[string]`` **DOMAIN-fractional_areas** Subgrid
      model fractional areas, see note above. `[-]`
    * `"snow death rate key`" ``[string]`` **DOMAIN-death_rate** Deals with last
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
  - `"surface-vapor_pressure_air`" `[Pa]`
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

Note that this PK is slaved to the flow PK.  This PK must be advanced first,
and be weakly coupled to the flow PK (or an MPC that advances the flow PK), and
the timestep of this PK must match that of the flow PK (e.g. do not try to
subcyle one or the other).  This uses a rather hacky, unconventional use of
time tags, where we use saturations and porosities at the NEXT time, but ASSUME
they are actually the values of the CURRENT time.  This saves having to stash a
copy of these variables at the CURRENT time, which would otherwise not be used.

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
  when it works, but it almost never works in real problems, except in
  column-based models, where it is perfect.

- "mstk implementation" MSTK implements an iterated, local optimization method
  that, one-at-a-time, moves nodes to try and match the volumes.  This has
  fewer issues with overfitting, but doesn't always do sane things, and can be
  expensive if iterations don't work well.  This is not particularly robust
  either, but it seems to be the preferred method for 2D/3D problems.

- "global optimization" attempts to directly form and solve the minimization
  problem to find the nodal changes that result in the target volumetric
  changes.  Note this has issues with overfitting, so penalty methods are used
  to smooth the solution of the problem.  This is currently disabled.

NOTE: all deformation options are treated EXPLICITLY, and depend only upon
values from the old time.

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

    - `"saturation_ice`" **DOMAIN-saturation_ice**
    - `"saturation_liquid`" **DOMAIN-saturation_liquid**
    - `"saturation_gas`" **DOMAIN-saturation_gas**
    - `"porosity`" **DOMAIN-porosity**
    - `"cell volume`" **DOMAIN-cell_volume**

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

   * `"evaluators`" ``[evaluator-typedinline-spec-list]`` A list of evaluators.

   * `"initial conditions`" ``[list]`` A list of constant-in-time data.  Note
     that `"initial conditions`" is not a particularly descriptive name here --
     PDE initial conditions are generally not here.  This list consists of

.. _evaluator-typedinline-spec:
.. admonition:: evaluator-typedinline-spec

   * `"evaluator type`" ``[string]`` Type of the evaluator Included for
     convenience in defining data that is not in the dependency graph,
     constants are things (like gravity, or atmospheric pressure) which are
     stored in state but never change.  Typically they're limited to scalars
     and dense, local vectors.


Example:

.. code-block:: xml

    <ParameterList name="state">
      <Parameter name="initialization filename" type="string" value="_CHECK00123.h5"/>
      <ParameterList name="evaluators">
        <ParameterList name="pressure">
          <Parameter name="evaluator type" type="string" value="primary variable" />
        </ParameterList>
      </ParameterList>

      <ParameterList name="initial conditions">
        <Parameter name="time" type="double" value="0.0">
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

 

Evaluators
==========

Evaluators are individual terms used to build up a PK or MPCs.  Each
term represents a variable in the equation, and can consist of primary
variables (those that are solved for by a PK solver), independent
variables (those that have no dependent variables but are provided by
the user as data), and secondary variables (those that are functions
of other variables).  Note that all three may be variable in space
and/or time.

.. contents:: **List of Evalutors**
   :local:

Primary Variables
-----------------


An evaluator with no dependencies that serves as the primary variable to be
solved for by a PK.  Note that users almost never are required to write an
input spec for these -- they are controlled by the PK and therefore the input
spec for this evaluator is written by that PK.

.. _evaluator-primary-spec:
.. admonition:: evaluator-primary-spec

   * `"tag`" ``[string]`` **""** Time tag at which this primary variable is used.




Independent Variables
---------------------

Independent variables are variables which do not depend upon other
dependent variables but are provided as data by the user.  Examples
include material properties, forcing datasets, etc.  These can be
provided in a few forms:

Constant
^^^^^^^^
 A field evaluator with no dependencies, a constant value.

This evaluator is typically used for providing data that is a simple constant
value.

This evaluator is used by providing the option:

`"evaluator type`" = `"independent variable constant`"

.. _independent-variable-constant-evaluator-spec:
.. admonition:: independent-variable-constant-evaluator-spec

   * `"value`" ``[double]`` The value.




From Function
^^^^^^^^^^^^^


This evaluator is typically used for providing data that are functions of space
and time.  The evaluator consists of a list of region,function pairs, and the
functions are evaluated across that region at each timestep.  If the problem is
time-independent, the `"constant in time`" option results in a performance
boost (as the functions need only be evaluated once).  This leverages the
exaustive functional format capability provided in Amanzi's Functions_ library.

This evaluator is used by providing the option:

`"evaluator type`" == `"independent variable`"

.. _independent-variable-function-evaluator-spec:
.. admonition:: independent-variable-function-evaluator-spec

   * `"constant in time`" ``[bool]`` **false** If true, only evaluate the
     functions once as they are time-independent.
   * `"function`" ``[composite-vector-function-spec-list]``




From File
^^^^^^^^^
 An evaluator with no dependencies specified by discrete data in a file.

This evaluator is typically used for providing data that are functions of space
and time.  Data is provided, discretely (e.g. with one data point per
cell/face/node), at a series of time slices.  The time slices are interpolated
linearly in time to provide the value.

Within the file, data is expected to meet the following (HDF5) layout::

   /time : a 1D array of length NTIMES, providing the time in seconds.
   /variable_name.ENTITY.DOF  (group)

      /0 : a 1D array of length NENTITIES, providing the values for each entity
           at time /time[0]
      /1 : ...
      /NTIMES-1 : 1D array at time /time[NTIMES-1]

This evaluator is used by providing the option:

`"evaluator type`" == `"independent variable`"

.. _independent-variable-from-file-evaluator-spec:
.. admonition:: independent-variable-from-file-evaluator-spec

   * `"filename`" ``[string]`` Path to the file.
   * `"variable name`" ``[string]`` Name of the dataset to read from the file.
   * `"domain name`" ``[string]`` **domain** Name of the domain on which the
      field is defined.
   * `"component name`" ``[string]`` **cell** Name of the component in the
     field to populate.
   * `"mesh entity`" ``[string]`` **cell** Name of the entity on which the
     component is defined.
   * `"number of dofs`" ``[int]`` **1** Number of degrees of freedom to read.
   * `"time function`" ``[function-spec]`` **optional** If provided, time is
     first manipulated by this function before interpolation.  This is useful
     for things like cyclic data, which can use a modulo time function to
     repeat the same data.

.. code-block:: xml

  <ParameterList name="field_evaluators">  <!-- parent list -->
  <ParameterList name="porosity">
    <Parameter name="field evaluator type" type="string" value="independent variable from file"/>
    <Parameter name="filename" type="string" value="_DATA_FILE.h5"/>
    <Parameter name="domain name" type="string" value="domain"/>
    <Parameter name="variable name" type="string" value="porosity"/>
    <Parameter name="component name" type="string" value="cell"/>
    <Parameter name="mesh entity" type="string" value="cell"/>
    <Parameter name="number of dofs" type="int" value="1"/>

    <ParameterList name="time function">  
      <Parameter name="times" type="Array(double)" value="{1.0, 2.0, 3.0}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

The field *porosity* is defined as a cell-based variable and
interpolated between three time intervals.





Secondary Variables
-------------------

All other evaluators are secondary variable evaluators, and these are
grouped by physics concept or process type.

Secondary variables, by definition, define functions that evaluate one
or more variables as a function of one or more variables.  Therefore
all secondary evaluators provide at least one "Key," which is the
variable(s) computed, and at least one "Dependency."

If the evaluator computes one and only one key, that key is provided
by the name of the parameter list in the `"evalutors`" list of State_.
If more than one Key is computed, then the second must either be
guessed by the code (for instance if the provided key is
"saturation_liquid", it is likely that the other key is
"saturation_gas") or provided by the user.  If more than one key is
computed, all of the keys computed can be specified exactly via the
input spec.  Keys are provided in one of two parameters:

* `"my variable key`" ``[string]`` Specifically name the variable used as "my variable"

* `"my variable key suffix`" ``[string]`` Name a suffix, and the
  variable is given by DOMAIN-SUFFIX, where the DOMAIN is given by the
  prefix in the evaluator list's name.  This is particularly useful
  for collections of enumerated PKs, e.g. columnar PKs, where the
  DOMAIN might be computed on the fly based on a column ID.

Dependencies use the same approach -- each dependency variable name
may include a default, and looks for a "key" and "key suffix" as
potential options.

As an example, a saturation evaluator may depend on pressure, and may
detail all of its names via something like:

.. code-block:: xml

    <ParameterList name="domain:1-saturation_liquid">
      <Parameter name="saturation gas key" type="string" value="domain:1-saturation_gas" />

      <Parameter name="pressure key suffix" type="string" value="pressure" />
      <!-- OR EQUIVALENTLY -->
      <Parameter name="pressure key" type="string" value="domain:1-pressure" />
    </ParameterList>


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
        <Parameter name="evaluator dependencies" type="Array(string)" value="{snow-cell_volume, snow-water_equivalent, surface-molar_density_liquid}" />
        <Parameter name="units" type="string" value="mol" />
      </ParameterList>

Multiplicative evaluator for `canopy-water_content`:

.. code-block:: xml

      <ParameterList name="canopy-water_content" type="ParameterList">
        <Parameter name="field evaluator type" type="string" value="multiplicative evaluator" />
        <Parameter name="evaluator dependencies" type="Array(string)" value="{canopy-cell_volume, canopy-water_equivalent, surface-molar_density_liquid}" />
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
 Capillary pressure for gas on a liquid.

.. _pc-liquid-evaluator-spec:
.. admonition:: pc-liquid-evaluator-spec

   KEYS:

   - `"pressure`" **DOMAIN-pressure**




Capillary pressure of liquid on ice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Water Retention Model and Relative Permeability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 Evaluates saturation through water retention models.

Water Retention Models (WRMs) determine the saturation as a function of
pressure and the relative permeability as a function of saturation.  Most
commonly used in practice is the van Genuchten model, but others are available default default;

`"evaluator type`" = `"wrm`"

.. _wrm-evaluator-spec:
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

.. _wrm-partition-typedinline-spec:
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

`"evaluator type`" = `"WRM rel perm`"

.. _rel-perm-evaluator-spec:
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
   - `"saturation_liquid`"
   - `"density`" (if `"use density on viscosity in rel perm`" == true)
   - `"viscosity`" (if `"use density on viscosity in rel perm`" == true)
   - `"surface relative permeability`" (if `"boundary rel perm strategy`" == `"surface rel perm`")




Van Genuchten Model
~~~~~~~~~~~~~~~~~~~
 WRMVanGenuchten : water retention model using van Genuchten's parameterization

van Genuchten's water retention curve.

.. _WRM-van-Genuchten-spec:
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

.. _wrm-implicit-permafrost-spec:
.. admonition:: wrm-implicit-permafrost-spec

    * `"converged tolerance`" ``[double]`` **1.e-12** Convergence tolerance of the implicit solve.
    * `"max iterations`" ``[int]`` **100** Maximum allowable iterations of the implicit solve.
    * `"solver algorithm [bisection/toms]`" ``[string]`` **bisection** Use bisection or the TOMS algorithm from boost.




Freezing point depression model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Freezing point depression, smoothed model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Painter's permafrost model with freezing point depression, smoothed.

.. _wrm-fpd-smoothed-permafrost-spec:
.. admonition:: wrm-fpd-smoothed-permafrost-spec

    * `"reference temperature [K]`" ``[double]`` **273.15** The phase transition point
    * `"interfacial tension ice-water [mN m^-1]`" ``[double]`` **33.1**
    * `"interfacial tension air-water [mN m^-1]`" ``[double]`` **72.7**
    * `"smoothing width [K]`" ``[double]`` **1.** Smoothing out the freeze curve allows this to be slightly easier to solve.
    * `"latent heat [J kg^-1]`" ``[double]`` **3.34e5** Latent heat of fusion
    * `"water density [kg m^-3]`" ``[double]`` **998.87** Density of water.  Note this probably should use the calculated value.

 


Interfrost model
~~~~~~~~~~~~~~~~


Sutra-ICE model
~~~~~~~~~~~~~~~


This model is based on the emperical freezing curve used by Sutra-Ice,
documented in papers by Voss & Walvoord.

.. _wrm-sutra-permafrost-model-spec:
.. admonition:: wrm-sutra-permafrost-model-spec

   * `"temperature transition [K]`" ``[double]`` thickness of the transition from frozen to thawed
   * `"residual saturation [-]`" ``[double]`` Standard residual saturation
   * `"freezing point [K]`" ``[double]`` **273.15**


 


Compressible porosity
^^^^^^^^^^^^^^^^^^^^^


Compressible grains are both physically realistic (based on bulk modulus) and a
simple way to provide a non-elliptic, diagonal term for helping solvers to
converge.

`"evaluator type`" = `"compressible porosity`"

.. _compressible-porosity-evaluator-spec:
.. admonition:: compressible-porosity-evaluator-spec

   * `"compressible porosity model parameters`" ``[compressible-porosity-standard-model-spec-list]``

   KEYS:

   - `"pressure`" **DOMAIN-pressure**
   - `"base porosity`" **DOMAIN-base_porosity**




Standard model
~~~~~~~~~~~~~~
 A simple model for allowing porosity to vary with pressure.

Based on a linear increase, i.e.

.. math::

   \phi = \phi_{base} + H(p - p_{atm}) * \alpha

where :math:`H` is the heaviside function and :math:`\alpha` is the provided
compressibility.  If the inflection point is set to zero, the above function is
exact.  However, then the porosity function is not smooth (has discontinuous
derivatives), so the inflection point smooths this with a quadratic that
matches the value and derivative at the inflection point and is 0 with 0 slope
at atmospheric pressure.

.. _compressible-porosity-standard-model-spec:
.. admonition:: compressible-porosity-standard-model-spec

   * `"region`" ``[string]`` Region on which this is applied.
   * `"pore compressibility [Pa^-1]`" ``[double]``  :math:`\alpha` as described above
   * `"pore compressibility inflection point [Pa]`" ``[double]`` **1000**

  The inflection point above which the function is linear.

Example:

.. code-block:: xml

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

.. _compressible-porosity-leijnse-model-spec:
.. admonition:: compressible-porosity-leijnse-model-spec

   * `"pore compressibility [Pa^-1]`" ``[double]`` :math:`\alpha` as described above
   * `"pore compressibility inflection point [Pa]`" ``[double]`` **1000** The
     inflection point above which the function is linear.

NOTE: additionally the user should provide a parameter in the `EWC
Globalization Delegate`_ to turn Leijnse model ON in the EWC calculations.

.. code-block:: xml

   <Parameter name="porosity leijnse model" type="bool" value="true"/>





Viscosity of water
^^^^^^^^^^^^^^^^^^

Two main viscosity models are commonly used -- a constant and one
which is temperature-dependent.  The viscosity of water is strongly
temperature dependent, so it is highly recommended to use that one if
the problem is nonisothermal.

Constant
~~~~~~~~
Like any quantity, a viscosity can simply be a constant value, at which
point it is not a secondary variable but an independent variable.

.. code-block:: xml

      <ParameterList name="viscosity_liquid" type="ParameterList">
        <Parameter name="field evaluator type" type="string" value="independent variable constant" />
        <Parameter name="value" type="double" value="8.9e-4" />
        <Parameter name="units" type="string" value="Pa s" />
      </ParameterList>

Nonisothermal
~~~~~~~~~~~~~


A non-isothermal viscosity model intended for use within a range of
temperatures from well below freezing to ~100C.

.. _viscosity-evaluator-spec:
.. admonition:: viscosity-evaluator-spec

   * `"viscosity model parameters`" ``[viscosity-typedinline-spec-list]``

   KEYS:

   - `"temperature`"

 



A constitutive relation for the viscosity of water as a function of temperature
in K, given as an empirical series expansion fit to data.

Used by setting

`"viscosity relation type`" = `"liquid water`"

.. _viscosity-water-spec:
.. admonition:: viscosity-water-spec

   NONE





Surface flow evaluators
-----------------------

Assorted evaluators used for surface flow, including potential
surfaces, Manning's conductivity, and their frozen equivalents.

Like the subsurface flow evaluators, many of these evaluators show up
in nearly all ATS simulations.  For real examples, see `ats-demos
<https://github.com/amanzi/ats-demos>`_


Ponded Depth or Water Height
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Computes ponded depth from surface water pressure.

.. math::
   h = \frac{H(p - p_{atm})}{\rho g}

where :math:`H` is the Heaviside function.

`"evaluator type`" = `"ponded depth`"

.. _height-evaluator-spec:
.. admonition:: height-evaluator-spec

   KEYS:

   - `"mass density`"
   - `"pressure`"




Ponded Depth, Frozen
^^^^^^^^^^^^^^^^^^^^



Effective, or Smoothed Height
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Computes ponded depth from surface water pressure using a smoothed term to make
derivative smooth near 0.  This is pretty much never used anymore.

`"evaluator type`" = `"effective height`"

.. _effective-height-evaluator-spec:
.. admonition:: effective-height-evaluator-spec

   * `"smoothing width [m]`" ``[double]`` **0.01** the width over which smoothing
     is applied.

   KEYS:

   - `"height`" The unsmoothed ponded depth




Unfrozen fraction
^^^^^^^^^^^^^^^^^


An empirical equation for freezing ponded water -- this is simply a smooth
sinusoidal curve from 0 to 1 over a given transition in temperature.

`"evaluator type`" = `"unfrozen fraction`"

.. _unfrozen-fraction-evaluator-spec:
.. admonition:: unfrozen-fraction-evaluator-spec

   * `"unfrozen fraction model`" ``[unfrozen-fraction-model-spec]``

   KEYS:

   - `"temperature`"




.. _unfrozen-fraction-model-spec:
.. admonition:: unfrozen-fraction-model-spec

   * `"transition width [K]`" ``[double]`` **0.2** Degrees over which to
     transition from no ice to all ice.

   * `"freezing point [K]`" ``[double]`` **273.15** Center of the transition,
     at this point unfrozen fraction is 0.5.

   * `"minimum unfrozen fraction [-]`` ``[double]`` **0** Sets a minimum value.




Unfrozen Flowing Depth
^^^^^^^^^^^^^^^^^^^^^^
 Evaluates the unfrozen mobile depth.

In freezing conditions, water is only mobile if it is unfrozen.  This evaluator
determines how much water is allowed to flow given that it is partially frozen.

.. math:

   \delta_{mobile} = \delta \chi^{\alpha}

Given a ponded depth, an unfrozen fraction, and an optional power-law exponent,
which we call the ice retardation exponent.

.. _unfrozen-effective-depth-evaluator-spec:
.. admonition:: unfrozen-effective-depth-evaluator-spec

  * `"ice retardation exponent [-]`" ``[double]`` **1.0** exponent alpha
    controlling how quickly ice turns off flow.

  DEPENDENCIES:
  - `"depth`" **DOMAIN-depth**
  - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**




SurfacePotential
^^^^^^^^^^^^^^^^^^^
 Evaluates the potential surface upon which overland flow acts.

.. math::
   h + z

`"evaluator type`" = 

.. _pres-elev-evaluator-spec:
.. admonition:: pres-elev-evaluator-spec

   KEYS:

   - `"height`" **DOMAIN-ponded_depth** Names the height variable. [m]
   - `"elevation`" **DOMAIN-elevation** Names the elevation variable. [m]


NOTE: This is a legacy evaluator, and is not in the factory, so need not be in
the input spec.  However, we include it here because this could easily be
abstracted for new potential surfaces, kinematic wave, etc, at which point it
would need to be added to the factory and the input spec.

NOTE: This could easily be replaced by a generic Additive_ Evaluator.




Overland Conductivity, sheet flow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. _`Overland Conductivity Evaluator`:
 Evaluates the conductivity of surface flow.

This implements the conductivity of overland flow, which is the nonlinear
coefficient in the diffusion wave equation.  The term is given by:

.. math:
   k = \frac{\delta^{coef}}{n_{mann} \sqrt(| \nabla z |)}

Optionally, this may include a density factor, typically a molar density, which
converts the flow law to water flux rather than volumetric flux.

Also, this evaluator can be used in snow redistribution, and in that case needs
some extra factors (timestep size) to ensure the correct flow law in that case.

`"evaluator type`" = `"overland conductivity`"

.. _overland-conductivity-evaluator-spec:
.. admonition:: overland-conductivity-evaluator-spec

   * `"include density`" ``[bool]`` **true** Include the density prefactor,
     converting the flux from volumetric flux to water flux.
   * `"dt factor [s]`" ``[double]`` **-1** The artificial timestep size used in calculating
      snow redistribution, only used in that case.
   * `"swe density factor [-]`" ``[double]`` **10** Ratio of water to snow density.

   DEPENDENCIES:

   - `"mobile depth`" **DOMAIN-mobile_depth** Depth of the mobile water; delta
     in the above equation.
   - `"slope`" **DOMAIN-slope_magnitude** Magnitude of the bed surface driving
     flow; | \nabla z | above.
   - `"coefficient`" **DOMAIN-manning_coefficient** Surface roughness/shape
     coefficient; n_{mann} above.
   - `"molar density liquid`" **DOMAIN-molar_density_liquid** If `"include
     density`" is true, the density.





Overland Conductivity, litter resistance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


(missing documentation!)






Thermodynamic evaluators
-------------------------

Internal energy, enthalpy, thermal conductivity, etc used for both
surface and subsurface transport of energy.

Internal energy
^^^^^^^^^^^^^^^


Computes (specific) internal energy of as a function of temperature.

`"evaluator type`" = `"iem`"

.. _iem-evaluator-spec:
.. admonition:: iem-evaluator-spec

   * `"IEM parameters`" ``[IEM-model-typedinline-spec-list]``

   KEYS:

   - `"temperature`"




Linear
~~~~~~
 Internal energy based on a linear fit.

Linear internal energy model -- function of Cv and temperature

.. math::

    u = L_f +  C_v * (T - T_{ref})

`"IEM type`" = `"linear`"

.. _IEM-model-linear-spec:
.. admonition:: IEM-model-linear-spec

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

`"IEM type`" = `"quadratic`"

.. _IEM-model-quadratic-spec:
.. admonition:: IEM-model-quadratic-spec

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


Computes (specific) internal energy of as a function of temperature and molar
fraction of water vapor in the gaseous phase.

`"evaluator type`" = `"iem water vapor`"

.. _iem-water-vapor-evaluator-spec:
.. admonition:: iem-water-vapor-evaluator-spec

   * `"IEM parameters`" ``[IEM-water-vapor-model-spec]``

   KEYS:

   - `"temperature`"
   - `"vapor molar fraction`"



 Internal energy model for air and water vapor.

.. _iem-water-vapor-model-spec:
.. admonition:: iem-water-vapor-model-spec

    * `"latent heat [J mol^-1]`" ``[double]`` Latent heat of vaporization
    * `"heat capacity [J mol^-1 K^-1]`" ``[double]`` C_v




Enthalpy
~~~~~~~~


Computes enthalpy [MJ / mol] of as a function of internal energy, pressure, and density.

.. math::
   e = u + 10^{-6} * \frac{p}{n_l}

`"evaluator type`" = `"enthalpy`"

.. _enthalpy-evaluator-spec:
.. admonition:: enthalpy-evaluator-spec

   * `"include work term`" ``[bool]`` **false** If false, e = u, ignoring the work term.

   KEYS:

   - `"internal energy`"
   - `"pressure`"
   - `"mass density`"




Thermal Conductivity, two phases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Thermal conductivity based on two-phases (air,liquid) composition of the
porous media.

`"evaluator type`" = `"two-phase thermal conductivity`"

.. _thermal-conductivity-twophase-evaluator-spec:
.. admonition:: thermal-conductivity-twophase-evaluator-spec

   * `"thermal conductivity parameters`" ``[thermal-conductivity-twophase-typedinline-spec-list]``

   KEYS:

   - `"porosity`"
   - `"saturation liquid`"




Wet-Dry Model
~~~~~~~~~~~~~


Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.

`"thermal conductivity type`" = `"two-phase wet/dry`"

.. _thermal-conductivity-twophase-wetdry-spec:
.. admonition:: thermal-conductivity-twophase-wetdry-spec

   * `"region`" ``[string]`` Region name on which to apply these parameters.
   * `"thermal conductivity, wet [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of saturated soil
   * `"thermal conductivity, dry [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of dry soil
   * `"unsaturated alpha [-]`" ``[double]`` Interpolating exponent
   * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

Example:

.. code:: xml

  <ParameterList name="thermal conductivity model">
    <Parameter name="thermal conductivity type" type="string" value="two-phase wet/dry"/>
    <Parameter name="thermal conductivity, wet [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity, dry [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="epsilon" type="double" value="1.e-10"/>
    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
  </ParameterList>

Units: [W m^-1 K^-1]




Peters-Lidard Model
~~~~~~~~~~~~~~~~~~~


A two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Emperical fit for dry conductivity from Peters-Lidard et al '98.

See Atchley et al GMD 2015 Supplementary Material for equations.

`"thermal conductivity type`" = `"two-phase Peters-Lidard`"

.. _thermal-conductivity-twophase-peterslidard-spec:
.. admonition:: thermal-conductivity-twophase-peterslidard-spec

    * `"region`" ``[string]`` Region name on which to apply these parameters.
    * `"thermal conductivity of soil [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of soil grains (not bulk soil)
    * `"thermal conductivity of liquid [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of liquid (water)
    * `"thermal conductivity of gas [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of gas (air)
    * `"unsaturated alpha [-]`" ``[double]`` Interpolating exponent
    * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

Example:

.. code:: xml

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="thermal conductivity type" type="string" value="two-phase Peters-Lidard"/>
    <Parameter name="thermal conductivity of soil [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of gas [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
    <Parameter name="epsilon" type="double" value="1.e-10"/>
  </ParameterList>

Units: [W m^-1 K^-1]




Thermal Conductivity, three phases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Thermal conductivity based on a three-phase (air,liquid,ice) composition of the
porous media.

`"evaluator type`" = `"three-phase thermal conductivity`"

.. _thermal-conductivity-threephase-evaluator-spec:
.. admonition:: thermal-conductivity-threephase-evaluator-spec

   * `"thermal conductivity parameters`" ``[thermal-conductivity-threephase-typedinline-spec-list]``

   KEYS:

   - `"porosity`"
   - `"saturation liquid`"
   - `"second saturation`"
   - `"temperature`"




Wet-Dry Model
~~~~~~~~~~~~~
 Three-phase thermal conductivity based on paper by Peters-Lidard.

A three-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Empirical relationship for frozen soil based on Peters-Lidard

See Atchley et al GMD 2015 Supplementary Material for equations.

`"thermal conductivity type`" = `"three-phase wet/dry`"

.. _thermal-conductivity-threephase-wetdry-spec:
.. admonition:: thermal-conductivity-threephase-wetdry-spec

   * `"region`" ``[string]`` Region name on which to apply these parameters.
   * `"thermal conductivity, saturated (unfrozen) [W m^-1 K^-1]`" ``[double]``
     Thermal conductivity of fully saturated, unfrozen bulk soil.
   * `"thermal conductivity, dry [W m^-1 K^-1]`" ``[double]`` Thermal conductivity
     of fully dried bulk soil.
   * `"unsaturated alpha unfrozen [-]`" ``[double]`` Interpolating exponent
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

`"thermal conductivity type`" = `"three-phase Peters-Lidard`"

.. _thermal-conductivity-threephase-peterslidard-spec:
.. admonition:: thermal-conductivity-threephase-peterslidard-spec

    * `"region`" ``[string]`` Region name on which to apply these parameters.
    * `"thermal conductivity of soil [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of soil grains (not bulk soil)
    * `"thermal conductivity of liquid [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of liquid (water)
    * `"thermal conductivity of gas [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of gas (air)
    * `"thermal conductivity of ice [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of ice
    * `"unsaturated alpha unfrozen [-]`" ``[double]`` Interpolating exponent
    * `"unsaturated alpha frozen [-]`" ``[double]`` Interpolating exponent
    * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

Example:

.. code:: xml

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

`"thermal conductivity type`" = `"three-phase volume averaged`"

See Atchley et al GMD 2015 Supplementary Material for equations.

.. _thermal-conductivity-volume-averaged-spec:
.. admonition:: thermal-conductivity-volume-averaged-spec

    * `"region`" ``[string]`` Region name on which to apply these parameters.
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


Thermal conductivity model with constant values as a function of temperature,
requires the sutra model for permafrost WRM to also be used.  This only exists
to support the INTERFROST comparison.

Usage:

.. code:: xml

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="Thermal Conductivity Type" type="string" value="sutra hacked"/>
    <Parameter name="thermal conductivity of frozen" type="double" value=""/>
    <Parameter name="thermal conductivity of mushy" type="double" value=""/>
    <Parameter name="thermal conductivity of unfrozen" type="double" value=""/>
    <Parameter name="residual saturation" type="double" value=""/>
  </ParameterList>

Units: [W m^-1 K^-1]



Thermal Conductivity, Surface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Thermal conductivity of surface water that can be either frozen or liquid phase.

`"evaluator type`" = `"surface thermal conductivity`"

.. _thermal-conductivity-surface-evaluator-spec:
.. admonition:: thermal-conductivity-surface-evaluator-spec

   * `"thermal conductivity parameters`" ``[thermal-conductivity-surface-spec]``

   KEYS:

   - `"unfrozen fraction`"
   - `"ponded depth`"

.. _thermal-conductivity-surface-spec:
.. admonition:: thermal-conductivity-surface-spec

   * `"thermal conductivity of water [W m^-1 K^-1]`" ``[double]`` **0.58**
   * `"thermal conductivity of ice [W m^-1 K^-1]`" ``[double]`` **2.18**
   * `"minimum thermal conductivity`" ``[double]`` **1.e-14**




Advected Energy Source
^^^^^^^^^^^^^^^^^^^^^^


Active Layer Averaged Temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Water Table
^^^^^^^^^^^


Computes the depth to a saturated water table.

`"evaluator type`" = `"water table depth`"

.. _water-table-depth-spec:
.. admonition:: water-table-depth-spec

    KEYS:

    - `"saturation of gas`" **SUBSURFACE_DOMAIN-saturation-gas**
    - `"subsurface cell volume`" **SUBSURFACE_DOMAIN-cell_volume**
    - `"surface cell volume`" **DOMAIN-cell_volume**





Thaw Depth
^^^^^^^^^^


Computes the depth to a saturated water table.

`"evaluator type`" = `"thaw depth`"

.. _thaw-depth-spec:
.. admonition:: thaw-depth-spec

    KEYS:

    - `"temperature`" **SUBSURFACE_DOMAIN-temperature**
    - `"subsurface cell volume`" **SUBSURFACE_DOMAIN-cell_volume**
    - `"surface cell volume`" **DOMAIN-cell_volume**





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

Standard Equation of State
^^^^^^^^^^^^^^^^^^^^^^^^^^
 EOSEvaluator is the interface between state/data and the model, an EOS.

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
"subgrid components," for instance snow covered and bare ground.  This
"subgrid" model allows smooth transitions between fully snow-covered
and fully bare ground.

These area fractions often get included as area-weights in calculating
full-cell quantities.

Two-component model
~~~~~~~~~~~~~~~~~~~~~~~~~~
 A subgrid model for determining the area fraction of snow vs not snow within a grid cell.

Uses a simple linear transition to vary between liquid and bare ground, and
another linear transition to vary between snow-covered and not-snow-covered.

Ordering of the area fractions calculated are: [bare ground/water, snow].

`"evaluator type`" = `"area fractions, two components`"

.. _area-fractions-twocomponent-evaluator-spec:
.. admonition:: area-fractions-twocomponent-evaluator-spec:

   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
         Mimimum area fraction allowed, less than this is rebalanced as zero.

   DEPENDENCIES:

   - `"snow depth`" ``[string]``

.. note:

   This evaluator also uses the LandCover_ types.  From that struct, it
   requires the value of the following parameters:

   - `"snow transition height [m]`" ``[double]`` **0.02**
      Minimum thickness for specifying the snow gradient.

.. note:

   This evaluator simplifies the situation by assuming constant density.  This
   make it so that ice and water see the same geometry per unit pressure, which
   isn't quite true thanks to density differences.  However, we hypothesize
   that these differences, on the surface (unlike in the subsurface) really
   don't affect the solution.




Three-component model
~~~~~~~~~~~~~~~~~~~~~~~~~~
 A subgrid model for determining the area fraction of land, water, and snow within a grid cell.

Uses a simple linear transition to vary between liquid and bare ground, and
another linear transition to vary between snow-covered and not-snow-covered.

Ordering of the area fractions calculated are: [bare ground, water, snow].

`"evaluator type`" = `"area fractions, three components`"

.. _area-fractions-threecomponent-evaluator-spec:
.. admonition:: area-fractions-threecomponent-evaluator-spec:

   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
      Mimimum area fraction allowed, less than this is rebalanced as zero.

   DEPENDENCIES:

   - `"snow depth`" **DOMAIN_SNOW-depth**
   - `"ponded depth`" **DOMAIN-ponded_depth**

.. note:

   This evaluator also uses the LandCover_ types.  From that struct, it
   requires the value of the following parameters:

   - `"snow transition height [m]`" ``[double]`` **0.02**
      Minimum thickness for specifying the snow gradient.
   - `"water transition height [m]`" ``[double]`` **0.02**
         Minimum thickness for specifying the water gradient.

.. note:

   This evaluator simplifies the situation by assuming constant density.  This
   make it so that ice and water see the same geometry per unit pressure, which
   isn't quite true thanks to density differences.  However, we hypothesize
   that these differences, on the surface (unlike in the subsurface) really
   don't affect the solution.




Three-component model, with microtopography
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 A subgrid model for determining the area fraction of land, open water, and snow within a grid cell.

Uses the subgrid equation from Jan et al WRR 2018 for volumetric or effective
ponded depth to determine the area of water, then heuristically places snow on
top of that surface.

`"evaluator type`" = `"area fractions, three components with microtopography`"

.. _area-fractions-threecomponent-microtopography-evaluator-spec:
.. admonition:: area-fractions-threecomponent-microtopography-evaluator-spec

   * `"snow transitional height [m]`" ``[double]`` **0.02**
     Minimum thickness for specifying the snow gradient.
   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
     Mimimum area fraction allowed, less than this is rebalanced as zero.
   * `"snow domain name`" ``[string]`` **DOMAIN_SNOW** A default is guessed at
     by replacing `"surface`" with `"snow`" in the this's domain.

   KEYS:

   - `"microtopographic relief`" **DOMAIN-microtopographic_relief**
     The name of del_max, the max microtopography value.
   - `"excluded volume`" **DOMAIN-excluded_volume**
     The name of del_excluded, the integral of the microtopography.
   - `"ponded depth`" **DOMAIN-pressure**
     The name of the surface water ponded depth.
   - `"snow depth`" **DOMAIN_SNOW-depth**
     The name of the snow depth.
   - `"volumetric snow depth`" **DOMAIN_SNOW-volumetric_depth**
     The name of the snow depth.


NOTE: this evaluator simplifies the situation by assuming constant density.
This make it so that ice and water see the same geometry per unit pressure,
which isn't quite true thanks to density differences.  However, we hypothesize
that these differences, on the surface (unlike in the subsurface) really don't
matter much. --etc





Potential Evapotranspiration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Models of potential evapotranspiration approximate the difference in
vapor pressure between the atmosphere and the soil as a function of
available energy, allowing the calculation of the max flux of ET that
the atmosphere could accept.  This can then be limited based on water
availability, etc.

Priestley-Taylor Potential Evapotranspiration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Evaluates potential evapotranpiration (PET) using Priestley & Taylor formulation.

This implementation is based on models provided in the PRMS-IV, Version 4, see
pages 90-93, Equations 1-57 to 1-60

Requires the following dependencies:

.. _pet-priestley-taylor-evaluator-spec:
.. admonition:: pet-priestley-taylor-evaluator-spec:

   * `"include limiter`" ``[bool]`` **false** If true, multiply potential ET by
     a limiter to get an actual ET.
   * `"limiter number of dofs`" ``[int]`` **1** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides how many.
   * `"limiter dof`" ``[int]`` **0** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides which one to use.
   * `"include 1 - limiter`" ``[bool]`` **false** If true, multiply potential
     ET by 1 - a limiter (e.g. a limiter that partitions between two pools) to
     get actual ET.
   * `"1 - limiter number of dofs`" ``[int]`` **1** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides how many.
   * `"1 - limiter dof`" ``[int]`` **0** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides which one to use.
   * `"sublimate snow`" ``[bool]`` **false** If true, use latent heat of
      vaporization of snow, not water.

   KEYS:

   - `"air temperature`" **DOMAIN-air_temperature** Air temp, in [K]
   - `"surface temperature`" **DOMAIN-temperature** Ground or leaf temp, in [K].  Note this may be the
      same as air temperature.
   - `"elevation`" **DOMAIN-elevation** Elevation [m]
   - `"net radiation`" **DOMAIN-net_radiation** [W m^-2] Net radiation balance, positive to the ground.
   - `"limiter`" [-] See `"include limiter`" above.
   - `"1 - limiter`" [-] See `"include 1 - limiter`" above.




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

Also known as Beta, or the water availability factor, or the plant wilting
factor, or the transpiration reduction function.

.. math::
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
   - `"saturation liquid`" **DOMAIN_SUB-saturation_liquid**
   - `"porosity`" **DOMAIN_SUB-porosity**
   - `"potential evaporation`" **DOMAIN_SUB-potential_evaporation**





Radiation Balance Terms
^^^^^^^^^^^^^^^^^^^^^^^

Often a balance of incoming and outgoing short and longwave radiations
are required to determine the energy available to go into latent heat,
and therefore potential evapotranspiration.

Note that incoming shortwave radiation is typically a user-provided
meterological forcing dataset.

Radiation Balance
~~~~~~~~~~~~~~~~~
 Evaluates a net radiation balance for ground and canopy.

Here the net radiation is positive for energy inputs to the layer.  Note that
ground is based on the two-channel (land + snow) while canopy is assumed to be
a simple, single layer.

Requires the use of LandCover types, for albedo and Beer's law coefficients.

This is combination of CLM v4.5 Tech Note and Beer's law for attenuation of
radiation absorption.  In particular, long-wave is exactly as Figure 4.1c in CLM
4.5 Tech Note.  The main difference comes in how absorptivity (which is equal
to emissivity, epsilon in that document) is defined.  Here we use Beer's law
which is an exponential decay with LAI.

Unlike CLM 4.5, here we do not split shortwave into direct and diffuse light.

Computes:

1. "surface radiation balance" -- Net radiation seen by the bare soil/ponded
   water, this includes radiation transmitted to the surface through the
   canopy, longwave emitted by the canopy, and less the longwave emitted by the
   surface itself.  [W m^-2] of actual area -- this does NOT include the
   surface area fraction factor which would be required to compute a total
   energy flux in W.
   
2. "snow radiation balance" -- Net radiation seen by the snow.  See surface
   above -- all are the same except using snow properties. [W m^-2]
   
3. "canopy radiation balance" -- this is a compute computation of the net
   radiation experienced by the canopy.  It includes the portion of shortwave
   and longwave from the atmosphere that are absorbed via Beer's law, minus the
   outgoing longwave emitted from the canopy, plus upward longwave radiation
   emitted by the snow and surface.  It also does not include any secondary
   bounces (e.g. reflected atmosphere->canopy->cloud->back to canopy, or
   transmitted by the canopy, reflected by snow/surface.

Requires the use of LandCover types, for canopy albedo and Beer's law
coefficients.

`"evaluator type`" = `"radiation balance, surface and canopy`"

.. _radiation-balance-evaluator-spec:
.. admonition:: radiation-balance-evaluator-spec

   KEYS:
   - `"surface albedos`" **SURFACE_DOMAIN-albedos**
   - `"surface emissivities`" **SURFACE_DOMAIN-emissivities**
   - `"incoming shortwave radiation`" **SURFACE_DOMAIN-incoming_shortwave_radiation**
   - `"incoming longwave radiation`" **SURFACE_DOMAIN-incoming_longwave_radiation**
   - `"surface temperature`" **SURFACE_DOMAIN-temperature**
   - `"snow temperature`" **SNOW_DOMAIN-temperature**
   - `"canopy temperature`" **CANOPY_DOMAIN-temperature**
   - `"leaf area index`" **CANOPY_DOMAIN-leaf_area_index**
   - `"area fractions`" **SURFACE_DOMAIN-area_fractions**

Note that this is a superset of the physics in the "canopy radiation
evaluator," and is therefore mutually exclusive with that model.
     



Canopy Radiation Balance
~~~~~~~~~~~~~~~~~~~~~~~~
 Evaluates the canopy radiation balance, providing canopy net and radiation to the snow/surface.

Computes and sums the downward radiation terms, determining the total radiation
sent down to the surface from the canopy and above.

Requires the use of LandCover types, for albedo and emissivity of the canopy
itself, along with Beer's law coefficients.

Computes:

1. canopy-downward_shortwave_radiation -- transmitted shortwave.  Note that
   incoming shortwave is attenuated by Beer's law, and partially transmitted
   without attenuation when there are gaps (e.g. LAI < 1) in the canopy.
   
2. canopy-downward_longwave_radiation -- transmitted longwave (see above,
   noting that Beer's law coefficients should be used that absorb most if not
   all the longwave radiation), along with longwave emitted by the canopy
   computed using a canopy leaf temperature and a Bolzmann equation.
   
3. canopy-downward_net_radiation -- this is a partial computation of the net
   radiation experienced by the canopy.  It includes the portion of shortwave
   and longwave from the atmosphere that are absorbed via Beer's law, minus the
   outgoing longwave emitted from the canopy (see downward above).  It does NOT
   include upward longwave radiation emitted by the snow or surface.  It also
   does not include any secondary bounces (e.g. reflected
   atmosphere->canopy->cloud->back to canopy, or transmitted by the canopy,
   reflected by snow).

Here the net radiation is positive for energy added to the canopy, while the
other two are positive for energy sent to the layer below.
   
In the canopy-downward_net_radiation, we cannot include the upward terms YET,
because these are a function of snow and surface temperature, which in turn
depend upon the downward radiation computed here.  So we choose to break the
loop here, by computing downard terms first, then iterating to compute snow
temperature, then compute upward terms.  The alternative would be to have a
formal snow energy PK that computed snow temperature, at which point we would
solve all of these balances to convergence simultaneously.

`"evaluator type`" = `"canopy radiation balance from above`"

.. _canopy-radiation-evaluator-spec:
.. admonition:: canopy-radiation-evaluator-spec

   KEYS:
   - `"incoming shortwave radiation`" **SURFACE_DOMAIN-incoming_shortwave_radiation**
   - `"incoming longwave radiation`" **SURFACE_DOMAIN-incoming_longwave_radiation**
   - `"canopy temperature`" **CANOPY_DOMAIN-temperature**
   - `"leaf area index`" **CANOPY_DOMAIN-leaf_area_index**

Note that this is a subset of the physics in the "radiation balance evaluator,"
and is therefore mutually exclusive with that model.




Surface Albedo
~~~~~~~~~~~~~~

Note that albedo is also a multiple subgrid component model, like
surface balance.

 Evaluates albedos and emissivities in a two-component subgrid model.

Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for two components -- snow and not snow
(water/ice/land).  Note this internally calculates albedo of snow based upon
snow density.

Components are indexed by: 0 = land/ice/water, 1 = snow.

Requires the use of LandCover types, for ground albedo and emissivity.

.. _albedo-evaluator-spec:
.. admonition:: albedo-evaluator-spec

   * `"albedo ice [-]`" ``[double]`` **0.44**
   * `"albedo water [-]`" ``[double]`` **0.1168**

   * `"emissivity ice [-]`" ``[double]`` **0.98**
   * `"emissivity water [-]`" ``[double]`` **0.995**
   * `"emissivity snow [-]`" ``[double]`` **0.98**

   KEYS:

   - `"subgrid albedos`" **DOMAIN-subgrid_albedos**
   - `"subgrid emissivities`" **DOMAIN-subgrid_emissivities**

   DEPENDENCIES:

   - `"snow density`" **SNOW_DOMAIN-density**
   - `"ponded depth`" **DOMAIN-ponded_depth**
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**




 Evaluates albedos and emissivities in a three-component subgrid model.

Evaluates the albedo and emissivity as an interpolation on the surface
properties and cover.  This allows for three components -- water/ice, land, and
snow.  Note this internally calculates albedo of snow based upon snow density.

Components are: 0 = land, 1 = ice/water, 2 = snow.

Requires the use of LandCover types, for ground albedo and emissivity.

.. _albedo-evaluator-subgrid-spec:
.. admonition:: albedo-evaluator-subgrid-spec

   * `"albedo ice [-]`" ``[double]`` **0.44**
   * `"albedo water [-]`" ``[double]`` **0.1168**

   * `"emissivity ice [-]`" ``[double]`` **0.98**
   * `"emissivity water [-]`" ``[double]`` **0.995**
   * `"emissivity snow [-]`" ``[double]`` **0.98**

   KEYS:

   - `"subgrid albedos`" **DOMAIN-subgrid_albedos**
   - `"subgrid emissivities`" **DOMAIN-subgrid_emissivities**

   DEPENDENCIES:

   - `"snow density`" **SNOW_DOMAIN-density**
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**




Incident Shortwave Radiation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Evaluates the radiation incident on a non-flat surface.

Aspect modified shortwave radiation is determined by a factor which is
multiplied by the 'incoming radiation incident on a flat surface' to determine
the 'incoming radiation incident on a sloping surface of a given aspect' as a
function of slope and aspect, Julian day of the year, and time of day.  The
latitude and Julian day of the year are used to modify this with both time of
day and seasonal changes of the planet.

Note that some careful checking and experimentation has found that, in
general, the daily average incoming radiation times the 12-noon aspect
modifier correlates reasonably well with the daily average of the
product of the hourly incoming radiation and the hourly aspect
modifier.  It is notably better than the daily average radiation times
the daily average aspect modifier.

This implementation is derived from `LandLab code
<https://github.com/landlab/landlab/blob/master/landlab/components/radiation/radiation.py>`_,
which is released under the MIT license.


.. _incident_shortwave_radiation_evaluator-spec:
.. admonition:: incident_shortwave_radiation_evaluator-spec

    * `"incident shortwave radiation parameters`" ``[incident_shortwave_radiation_model-spec]``

    KEYS:
    - `"slope`" **DOMAIN-slope_magnitude**
    - `"aspect`" **DOMAIN-aspect**
    - `"incoming shortwave radiation`" **DOMAIN-incoming_shortwave_radiation**



 Evaluates shortwave as a function of slope/aspect/etc.

.. _incident_shortwave_radiation_model-spec:
.. admonition:: incident_shortwave_radiation_model-spec

   * `"latitude [degrees]`" ``[double]`` Latitude of the site.  A single
     typical value is fine for most domains, even relatively large ones
     (e.g. HUC4).
   * `"daily averaged`" ``[bool]`` **true** Compute a daily averaged radiation, as
     opposed to a time-specific, diurnal-cycle-resolving value.
   * `"day of year at time 0 [Julian days]`" ``[int]`` **0** (Jan 1).  ATS has
     no notion of absolute time, so to do things that depend upon planetary
     dynamics we must know what the day of the year is.  Typically this is set
     by your meteorological data -- set this to be equal to the day of year of
     met data's time 0.




Longwave Radiation
~~~~~~~~~~~~~~~~~~
 Evaluates incoming longwave radiation from vapor pressure and air temperature.

.. _longwave_evaluator-spec:
.. admonition:: longwave_evaluator-spec


    DEPENDENCIES:

    * `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
    * `"vapor pressure air key`" ``[string]`` **DOMAIN-vapor_pressure_air**




Full Surface Energy Balance Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Finally, in addition to the potential-based models above, a few
full-physics model are available.  These are often based on older,
monolithic models.

Bare Soil Surface Energy Balance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Calculates source terms for surface fluxes to and from the atmosphere and a ground surface.

The ground is assumed to consist of two potential area-fraction components --
snow and no-snow.  In the case of snow on the ground, this solves for a snow
temperature, given a skin temperature, that satisfies a energy balance
equation.  In the case of no-snow, this calculates a conductive heat flux to
the ground from the atmosphere.

`"evaluator type`" = `"surface energy balance, two components`"

.. _seb-twocomponent-evaluator-spec:
.. admonition:: seb-twocomponent-evaluator-spec

   * `"wind speed reference height [m]`" ``[double]`` **2.0** Reference height at which
     wind speed is measured.
   * `"minimum wind speed [m s^-1]`" ``[double]`` **1.0** Sets a floor on wind speed for
     potential wierd data.  Models have trouble with no wind.

   * `"save diagnostic data`" ``[bool]`` **false** Saves a suite of diagnostic variables to vis.

   * `"surface domain name`" ``[string]`` **DEFAULT** Default set by parameterlist name.
   * `"subsurface domain name`" ``[string]`` **DEFAULT** Default set relative to surface domain name.
   * `"snow domain name`" ``[string]`` **DEFAULT** Default set relative to surface domain name.

   KEYS:

   - `"surface water source`" **DOMAIN-water_source**  [m s^-1]
   - `"surface energy source`" **DOMAIN-total_energy_source** [MW m^-2]
   - `"subsurface water source`" **DOMAIN-water_source**  [mol s^-1]
   - `"subsurface energy source`" **DOMAIN-total_energy_source** [MW m^-3]
   - `"snow mass source - sink`" **DOMAIN-source_sink** [m_SWE s^-1]
   - `"new snow source`" **DOMAIN-source** [m_SWE s^-1]

   - `"albedo`" **DOMAIN-albedo** [-] A single variate diagnostic of the final albedo.
   - `"snowmelt`" **DOMAIN_SNOW-melt** [m_SWE s^-1]
   - `"evaporation`" **DOMAIN-evaporative_flux** [m s^-1]
   - `"snow temperature`" **DOMAIN_SNOW-temperature** [K]
   - `"sensible heat flux`" **DOMAIN-qE_sensible_heat** [W m^-2]
   - `"latent heat of evaporation`" **DOMAIN-qE_latent_heat** [W m^-2]
   - `"latent heat of snowmelt`" **DOMAIN-qE_snowmelt** [W m^-2]
   - `"outgoing longwave radiation`" **DOMAIN-qE_lw_out** [W m^-2]
   - `"conducted energy flux`" **DOMAIN-qE_conducted** [W m^-2]

   DEPENDENCIES:

   - `"incoming shortwave radiation`" **DOMAIN-incoming_shortwave_radiation** [W m^-2]
   - `"incoming longwave radiation`" **DOMAIN-incoming_longwave_radiation** [W m^-2]
   - `"air temperature`" **DOMAIN-air_temperature** [K]
   - `"vapor pressure air`" **DOMAIN-vapor_pressure_air** [Pa]
   - `"wind speed`" **DOMAIN-wind_speed** [m s^-1]
   - `"precipitation rain`" **DOMAIN-precipitation_rain** [m s^-1]
   - `"precipitation snow`" **DOMAIN_SNOW-precipitation** [m_SWE s^-1]

   - `"snow depth`" **DOMAIN_SNOW-depth** [m]
   - `"snow density`" **DOMAIN_SNOW-density** [kg m^-3]
   - `"snow death rate`" **DOMAIN_SNOW-death_rate** [m s^-1]  Snow "death" refers to the last bit of snowmelt that we want to remove discretely.
   - `"ponded depth`" **DOMAIN-ponded_depth** [m]
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction** [-]  1 --> all surface water, 0 --> all surface ice
   - `"subgrid albedos`" **DOMAIN-albedos** [-] Dimension 2 field of (no-snow, snow) albedos.
   - `"subgrid emissivity`" **DOMAIN-emissivities** [-] Dimension 2 field of (no-snow, snow) emissivities.
   - `"area fractions`" **DOMAIN-fractional_areas** Dimension 2 field of (no-snow, snow) area fractions (sum to 1).

   - `"temperature`" **DOMAIN-temperature**  [K] surface skin temperature.
   - `"pressure`" **DOMAIN-pressure** [Pa] surface skin pressure.
   - `"gas saturation`" **DOMAIN_SS-saturation_gas** [-] subsurface gas saturation
   - `"liquid saturation`" **DOMAIN_SS-saturation_liquid** [-] subsurface liquid saturation
   - `"porosity`" [-] subsurface porosity
   - `"subsurface pressure`" **DOMAIN_SS-pressure** [Pa]
   - `"molar density liquid`" **DOMAIN-molar_density_liquid** [mol m^-3]
   - `"mass density liquid`" **DOMAIN-mass_density_liquid** [kg m^-3]


.. note:

   This also depends upon multiple parameters from the LandCover_ types:

   - `"roughness length of bare ground [m]`" ``[double]`` **0.04** Defines a fetch controlling
     latent and sensible heat fluxes.
   - `"roughness length of snow-covered ground [m]`" ``[double]`` **0.004** Defines a
     fetch controlling latent and sensible heat fluxes.
   - `"dessicated zone thickness [m]`" ``[double]`` Thickness of the immediate surface
     layer over which vapor pressure diffusion must move water to evaporate
     from dry soil.  More implies less evaporation.
   - `"snow transition depth [m]`" **0.02** Snow height at which bare
     ground starts to stick out due to subgrid topography, vegetation, etc.
     Defines a transitional zone between "snow-covered" and "bare ground".
   - `"water transition depth [m]`" **0.02** Ponded depth at which bare
     ground starts to stick out due to subgrid topography, vegetation, etc.
     Defines a transitional zone between "water-covered" and "bare ground".




 Calculates source terms for surface fluxes to and from the atmosphere and a ground surface.

Calculates source terms for surface fluxes to and from the atmosphere and a
ground surface characterized by three components -- snow, water-covered ground,
and vegetated/bare soil.

The surface energy balance on these area weighted patches are individually
calculated then averaged to form the total quantities.  All down- and
up-scaling of relevant quantities are done through the area weighting, which is
calculated by a minimum threshold in snow and a depression depth/geometry-based
approach for water.  All snow is assumed to first cover water (likely ice),
then cover land, as both water and snow prefer low-lying depressions due to
gravity- and wind-driven redistributions, respectively.

`"evaluator type`" = `"surface energy balance, two components`"

.. _seb-threecomponent-evaluator-spec:
.. admonition:: seb-threecomponent-evaluator-spec

   * `"wind speed reference height [m]`" ``[double]`` **2.0** Reference height at which
     wind speed is measured.
   * `"minimum wind speed [m s^-1]`" ``[double]`` **1.0** Sets a floor on wind speed for
     potential wierd data.  Models have trouble with no wind.

   * `"save diagnostic data`" ``[bool]`` **false** Saves a suite of diagnostic variables to vis.

   * `"surface domain name`" ``[string]`` **DEFAULT** Default set by parameterlist name.
   * `"subsurface domain name`" ``[string]`` **DEFAULT** Default set relative to surface domain name.
   * `"snow domain name`" ``[string]`` **DEFAULT** Default set relative to surface domain name.

   KEYS:

   - `"surface water source`" **DOMAIN-water_source**  [m s^-1]
   - `"surface energy source`" **DOMAIN-total_energy_source** [MW m^-2]
   - `"subsurface water source`" **DOMAIN-water_source**  [mol s^-1]
   - `"subsurface energy source`" **DOMAIN-total_energy_source** [MW m^-3]
   - `"snow mass source - sink`" **DOMAIN-source_sink** [m_SWE s^-1]
   - `"new snow source`" **DOMAIN-source** [m_SWE s^-1]

   - `"albedo`" **DOMAIN-albedo** [-] A single variate diagnostic of the final albedo.
   - `"snowmelt`" **DOMAIN_SNOW-melt** [m_SWE s^-1]
   - `"evaporation`" **DOMAIN-evaporative_flux** [m s^-1]
   - `"snow temperature`" **DOMAIN_SNOW-temperature** [K]
   - `"sensible heat flux`" **DOMAIN-qE_sensible_heat** [W m^-2]
   - `"latent heat of evaporation`" **DOMAIN-qE_latent_heat** [W m^-2]
   - `"latent heat of snowmelt`" **DOMAIN-qE_snowmelt** [W m^-2]
   - `"outgoing longwave radiation`" **DOMAIN-qE_lw_out** [W m^-2]
   - `"conducted energy flux`" **DOMAIN-qE_conducted** [W m^-2]

   DEPENDENCIES:

   - `"incoming shortwave radiation`" **DOMAIN-incoming_shortwave_radiation** [W m^-2]
   - `"incoming longwave radiation`" **DOMAIN-incoming_longwave_radiation** [W m^-2]
   - `"air temperature`" **DOMAIN-air_temperature** [K]
   - `"vapor pressure air`" **DOMAIN-vapor_pressure_air** [Pa]
   - `"wind speed`" **DOMAIN-wind_speed** [m s^-1]
   - `"precipitation rain`" **DOMAIN-precipitation_rain** [m s^-1]
   - `"precipitation snow`" **DOMAIN_SNOW-precipitation** [m_SWE s^-1]

   - `"snow depth`" **DOMAIN_SNOW-depth** [m]
   - `"snow density`" **DOMAIN_SNOW-density** [kg m^-3]
   - `"snow death rate`" **DOMAIN_SNOW-death_rate** [m s^-1]  Snow "death" refers to the last bit of snowmelt that we want to remove discretely.
   - `"ponded depth`" **DOMAIN-ponded_depth** [m]
   - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction** [-]  1 --> all surface water, 0 --> all surface ice
   - `"subgrid albedos`" **DOMAIN-albedos** [-] Dimension 2 field of (no-snow, snow) albedos.
   - `"subgrid emissivity`" **DOMAIN-emissivities** [-] Dimension 2 field of (no-snow, snow) emissivities.
   - `"area fractions`" **DOMAIN-fractional_areas** Dimension 2 field of (no-snow, snow) area fractions (sum to 1).

   - `"temperature`" **DOMAIN-temperature**  [K] surface skin temperature.
   - `"pressure`" **DOMAIN-pressure** [Pa] surface skin pressure.
   - `"gas saturation`" **DOMAIN_SS-saturation_gas** [-] subsurface gas saturation
   - `"liquid saturation`" **DOMAIN_SS-saturation_liquid** [-] subsurface liquid saturation
   - `"porosity`" [-] subsurface porosity
   - `"subsurface pressure`" **DOMAIN_SS-pressure** [Pa]
   - `"molar density liquid`" **DOMAIN-molar_density_liquid** [mol m^-3]
   - `"mass density liquid`" **DOMAIN-mass_density_liquid** [kg m^-3]


.. note:

   This also depends upon multiple parameters from the LandCover_ types:

   - `"roughness length of bare ground [m]`" ``[double]`` **0.04** Defines a fetch controlling
     latent and sensible heat fluxes.
   - `"roughness length of snow-covered ground [m]`" ``[double]`` **0.004** Defines a
     fetch controlling latent and sensible heat fluxes.
   - `"dessicated zone thickness [m]`" ``[double]`` Thickness of the immediate surface
     layer over which vapor pressure diffusion must move water to evaporate
     from dry soil.  More implies less evaporation.




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
 Computes the potential surface upon which snow redistribution acts.

Snow distribution moves snow precipitation around to ensure that low-lying
areas fill up first.  This represents a potential field driving movement of
snow precip, so defines what we mean by low-lying.

.. math::
    h + z + h_{{snow}} + dt * P_{{snow}}

`"evaluator type`" = `"snow skin potential`"

.. _snow-skin-potential-evaluator-spec:
.. admonition:: snow-skin-potential-evaluator-spec

   * `"dt factor [s]`" ``[double]`` A free-parameter factor for providing a
     time scale for diffusion of snow precipitation into low-lying areas.
     Typically on the order of 1e4-1e7. This timestep times the wave speed of
     snow provides an approximate length of how far snow precip can travel.
     Extremely tunable! [s]
                
   KEYS:
                
   - `"ponded depth`" **SURFACE_DOMAIN-ponded_depth** [m]
   - `"snow depth`" **DOMAIN-depth** [m]
   - `"precipitation snow`" **DOMAIN-precipitation** [m s^-1]
   - `"elevation`" **SURFACE_DOMAIN-elevation** [m]


Example:

.. code-block:: xml

  <ParameterList name="snow_skin_potential" type="ParameterList">
    <Parameter name="evaluator type" type="string" value="snow skin potential" />
    <Parameter name="dt factor [s]" type="double" value="864000.0" />
  </ParameterList>




Snow Melt Rate
^^^^^^^^^^^^^^
 Evaluates snow melt via USDA - Natural Resources Conservation Service model

From:  National Engineering Handbook (NEH) part 630, Chapter 11

Snow melt rate is given by:

.. math::
   SM = H(T_{snow}^{expected} - 273.15) R

where :math:`R` is the snow melt rate per degree-day and
:math:`T_{snow}^{expected}` is the expected snow temperature, which is
typically given by :math:`T_{snow}^{expected} = T_{air} - \Delta`, where
:math:`\Delta` is the expected air-snow temperature difference [C] (note this
is NOT prescribed here -- the user must supply the expected snow temperature
via an evalutor).

Note that the Heaviside function is used to ensure this is only active when the
expected snow temperature is above 0 C.

Then, a linear transition factor is applied to ensure that this goes to zero as
the snow SWE goes to zero.  That factor is 1 at the snow transition depth, and
0 when snow SWE is 0.  This uses LandCover for the snow_ground_transition
parameter.

.. _snow-meltrate-evaluator-spec:
.. admonition:: snow-meltrate-evaluator-spec

   * `"snow melt rate [mm day^-1 C^-1]`" ``[double]`` **2.74**
     the melt rate per degree-day, above 0 C, e.g. :math:`R` above.

   KEYS:

   - `"snow water equivalent`" **DOMAIN-water_equivalent**
   - `"expected snow temperature`"  **DOMAIN-expected_temperature**





Canopy evaluators
-----------------

Interception
^^^^^^^^^^^^
 Partitions water between interception and throughfall, combining throughfall with drainage.

Based on CLM 4.5 and Lawrence et al 2007, interception is given by:

.. math::
   I = (P_{rain} + P_{snow}) * \alpha * (1 - exp(-.5(LAI)))

Throughfall is given by:

.. math::
   T = (P_{rain} + P_{snow}) - I

Drainage is provided as input here, as a total drainage from the canopy.  The
phase of this drainage is assumed to match the phase of the precipitation.  So
if it is raining, drainage is rain, while if it is 50/50 rain and snow,
drainage is also 50/50 rain and snow.  If total precipitation is 0, then
drainage is partitioned by air temperature (above 0C --> all rain, otherwise
all snow).  This evaluator partitions the drainage and sums it with throughfall
to compute the total source, in each phase, to the layer below the canopy (snow
and/or ground surface).

.. _interception-fraction-evaluator-spec:
.. admonition:: interception-fraction-evaluator-spec

   * `"interception fraction parameters`" ``[interception-fraction-model-spec]``

   MY KEYS:
   - `"interception`" **DOMAIN-interception**
   - `"throughfall and drainage rain`" **DOMAIN-throughfall_drainage_rain**
   - `"throughfall and drainage snow`" **DOMAIN-throughfall_drainage_snow**

   KEYS:
   - `"area index`" **DOMAIN-area_index**
   - `"precipitation rain`" **DOMAIN_SURFACE-precipitation_rain**
   - `"precipitation snow`" **DOMAIN_SNOW-precipitation**
   - `"drainage`" **DOMAIN-drainage**
   - `"air temperature`" **DOMAIN_SURFACE-air_temperature**



 Fraction of incoming water that is intercepted.

Based on CLM 4.5 and Lawrence et al 2007:

.. math::
  f_I = \alpha * (1 - exp(-.5 LAI)))

The interception fraction is everything here after the precip.

.. _interception-fraction-model-spec:
.. admonition:: interception-fraction-model-spec

   * `"leaf area interception fraction [-]`" ``[double]`` **0.25** The alpha
      term, this describes the fraction of leaf area that intercepts water.




Drainage
^^^^^^^^
 Drainage rate from the canopy to the lower layers.

A simple model based on relaxation from current water content to a saturated water content.

.. code::
   
          |
          | source
          V
         /   \
      I /     \
       V       |
   --Theta--    | T
       ^       |
       | D     |
       V       V
   -- -- -- -- -- -- --


This is the model for drainage D.

Drainage is given by:

.. math::
   D = max(0, \frac{(\Theta - \Theta_sat)}{\tau})

.. _drainage-evaluator-spec:
.. admonition:: drainage-evaluator-spec

   * `"drainage timescale [s]`" ``[double]`` **864** Timescale over which drainage occurs.
   * `"saturated specific water content [m^3 H2O / m^2 leaf area]`" ``[double]`` **1e-4**
      The thickness of the wetting surface -- determines maximum canopy water storage.\

   KEYS:
   - "area index"
   - "water equivalent"







Carbon and Biogeochemistry Models
---------------------------------

Carbon Decomposition Rate
^^^^^^^^^^^^^^^^^^^^^^^^^


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

Note that this evaluator fills exactly one subdomain in a domain set -- there
will be N evaluators each filling one subdomain.

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
 Evaluates the elevation (z-coordinate) and slope magnitude of a mesh.

`"evaluator type`" = `"meshed elevation`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

.. _meshed-elevation-evaluator-spec:
.. admonition:: meshed-elevation-evaluator-spec

   * `"parent domain name`" ``[string]`` **SUBSURFACE_DOMAIN** Domain name of the parent
     mesh, which is the 3D version of this domain.  Attempts to generate an
     intelligent default by stripping "surface" from this domain.
   * `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the
     elevation changes in time, and adds the `"deformation`" dependency.
   
   MY KEYS:

   - `"elevation`" **DOMAIN-elevation** Name the elevation variable. [m]
   - `"slope magnitude`" **DOMAIN-slope_magnitude** Name the elevation variable. [-]

   KEYS:

   - `"deformation`" **optional** If `"dynamic mesh`" == True, this is required
     to tell when the mesh has been deformed.




Elevation of a Column
^^^^^^^^^^^^^^^^^^^^^
 ElevationEvaluatorColumn: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.

Evaluates elevation, slope, and aspect of the "surface_star" mesh of the Arctic
Intermediate Scale Model (ISM).

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h
z|`` Note that this is evaluator is different from both the standard elevation
evaluator, which would take elevation/slope/aspect from the parent 3D mesh, and
from the standard "single cell" column mesh elevation, which can do the same.
Instead, it is a mix of:

- elevation is aggregated from the column meshes' top faces

- slope is calculated using a cell-neighbor differencing scheme.  Unlike 3D
  meshes, in the ISM, we deform columns of cells and so do not have faces.
  Therefore, we use all neighboring cells to compute an average normal for the
  cell, and use that to compute slope and aspect.

- aspect is set to 0.  It could easily be calculated using the same normal as
  the slope algorithm, but is not done currently.

`"evaluator type`" = `"elevation column`"

.. _column-elevation-evaluator-spec:
.. admonition:: column-elevation-evaluator-spec

   * `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
   * `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation
     variable. [-]
   * `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the
     elevation changes in time, and adds the `"deformation`" and
     `"base_porosity`" dependencies.
   * `"parent domain name`" ``[string]`` **DOMAIN** Domain name of the parent
     mesh, which is the 3D version of this domain.  In the columnar meshes the
     surface elevation and slope are assigned based on the columns and not the
     base 3D domain.

Example:

.. code-block:: xml

  <ParameterList name="surface_star-elevation">
    <Parameter name="evaluator type" type="string" value="column elevation"/>
  </ParameterList>




Depth Evaluator
^^^^^^^^^^^^^^^
 Computes depth, positive downward relative to the surface, of mesh cells.

Computes cell depths only.

Note that two algorithms for computing the depth are available:

`"cell centroid`" uses the cell centroid, as computed by the Mesh, directly.

`"mean face centroid`" is the default, as the cell centroid can have problems
in the case of non-planar faces.  Instead, it uses the mean of the above and
below face centroids in place of the face centroid, with the implicit
assumption that dz is uniform (e.g. this is an extruded mesh).

`"evaluator type`" = `"depth`"

.. _depth-evaluator-spec:
.. admonition:: depth-evaluator-spec

   * `"algorithm`" ``[string]`` **mean face centroid** Valid is `"mean face centroid`"
     and `"cell centroid`", see above.

 


Cell Volume evaluator
^^^^^^^^^^^^^^^^^^^^^


Deforming Cell Volume evaluator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



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

   * `"enforce positivity`" ``[bool]`` **false** If true, max the result with 0.

   * `"DEPENDENCY coefficient`" ``[double]`` **1.0** A multiple for each dependency in
     the list below.

   ONE OF

   * `"dependencies`" ``[Array(string)]`` The things to sum.

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

   * `"DEPENDENCY dof`" ``[double]`` **0** Degree of Freedom for each given
     dependency to use in the multiplication.  NOTE, this should only be
     provided if the dependency has more than 1 DoF -- if it just has one
     leaving this blank results in better error checking than providing the
     value 0 manually.

   ONE OF

   * `"dependencies`" ``[Array(string)]`` The fields to multiply.

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
     variable's domain (typically "surface" or "surface_column:\*") and is
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

Uses functions to evaluate arbitrary algebraic functions of its dependencies.

For example, one might write a dependency:

  VARNAME = 0.2 * DEP1 - DEP2 + 3

`"evaluator type`" = `"secondary variable from function`"

.. _secondary-variable-from-function-evaluator-spec:
.. admonition:: secondary-variable-from-function-evaluator-spec

   ONE OF:

   * `"functions`" ``[composite-vector-function-spec-list]`` Note this is used
     for multiple Degress of Freedom.

   OR:

   * `"function`" ``[composite-vector-function-spec]`` Used for a single degree
     of freedom.

Example:

.. code:: xml

   <ParameterList name="VARNAME">
     <Parameter name="field evaluator type" type="string" value="algebraic variable from function"/>
     <Parameter name="evaluator dependencies" type="Array{string}" value="{DEP1, DEP2}"/>
     <ParameterList name="function">
       <ParameterList name="function-linear">
         <Parameter name="x0" type="Array(double)" value="{0.0,0.0}" />
         <Parameter name="y0" type="double" value="3." />
         <Parameter name="gradient" type="Array(double)" value="{0.2, -1}" />
       </ParameterList>
     </ParameterList>
   </ParameterList>






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
available options for initializing various data in state differ by
data type.



.. _constants-scalar-spec:
.. admonition:: constants-scalar-spec

   * `"value`" ``[double]`` Value of a scalar constant

.. _constants-dense-vector-spec:
.. admonition:: constants-dense-vector-spec

   * `"value`" ``[Array(double)]`` Value of a dense, local vector.

.. _constants-point-spec:
.. admonition:: constants-point-spec

   * `"value`" ``[Array(double)]`` Array containing the values of the point.


.. _constants-composite-vector-spec:
.. admonition:: constants-composite-vector-spec

   * `"constant`" ``[double]`` **optional** Constant value.
   * `"value`" ``[double]`` **optional** Constant value, same as `"constant`" above.
   * `"function`" ``[composite-vector-function-spec-list]`` **optional**
     Initialize from a function, see CompositeVectorFunction_
   * `"restart file`" ``[string]`` **optional** Path to a checkpoint file from
     which to read the values.
   * `"cells from file`" ``[string]`` **optional** Same as `"restart file`",
     but only reads the cell component.
   * `"exodus file initialization`" ``[exodus-file-initialization-spec]``
     **optional** See `Exodus File Initialization`_.
   * `"initialize from 1D column`" ``[column-file-initialization-spec]``
     **optional** See `Column File Initialization`_.




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
   <ParameterList name="water flux">
     <ParameterList name="BC north">
       <Parameter name="regions" type="Array(string)" value="{north}"/>
       <ParameterList name="outward water flux">
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


Neumann (water flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for both surface and subsurface flows, this provides water flux data (in [mol m^-2 s^-1], for the subsurface, or [mol m^-1 s^-1] for the surface, in the outward normal direction) on boundaries.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="water flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward water flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="-1.e-3"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Neumann (fix level flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for surface only,  this provides fixed level ([m])  velocity data (in [m s^-1], in the outward normal direction on boundaries.

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
                 <Parameter name="value" type="double" value="2.5"/>
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
       <ParameterList name="outward water flux">
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
       <Parameter name="bc types" type="Array(string)" value="{head, water flux}"/>
       <Parameter name="bc functions" type="Array(string)" value="{boundary head, outward water flux}"/>

       <ParameterList name="water flux">
         <ParameterList name="BC west">
           <Parameter name="regions" type="Array(string)" value="{surface west}"/>
           <ParameterList name="outward water flux">
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



Energy boundary conditions must follow the general format shown in
BoundaryConditions_.  Energy-specific conditions implemented include:

Dirichlet (temperature) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Provide temperature data in units of [K].

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="temperature">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary temperature">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="276.15."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Neumann (diffusive energy flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that for an energy equation, there are two mechanisms by which energy can
be fluxed into the domain: through diffusion and through advection.  This
boundary condition sets the diffusive flux of energy, and allows the advective
flux to be whatever it may be.  Frequently this is used in combination with
boundaries where water is expected to be advected out of the domain, and we
wish to allow the energy of that water to be advected away with it, but wish to
independently specify diffusive fluxes.  This can also be used in cases where
the mass flux is prescribed to be zero (e.g. bottom boundaries, where this
might be the geothermal gradient).

Units are in **[MW m^-2]**, noting the deviation from SI units!

Example:

.. code-block:: xml

  <ParameterList name="boundary conditions">
   <ParameterList name="diffusive flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward diffusive flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Note that another commonly implemented boundary condition is one where the
diffusive flux is prescribed, and also the temperature of incoming water is
prescribed.  This is not currently implemented, but would be straightforward to
do so if requested.


Neumann (total energy flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This boundary condition sets the total flux of energy, from both advection and
diffusion.  This is not used all that often in real applications, but is common
for benchmarks or other testing.

Units are in **[MW m^-2]**, noting the deviation from SI units!

Example:

.. code-block:: xml

  <ParameterList name="boundary conditions">
   <ParameterList name="enthalpy flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward enthalpy flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>











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

.. code-block:: xml

  <ParameterList name="time integrator">
    <Parameter name="time integration method" type="string" value="BDF1"/>
    <ParameterList name="BDF1">
      <Parameter name="max preconditioner lag iterations" type="int" value="5"/>
      <Parameter name="freeze preconditioner" type="bool" value="false"/>
      <Parameter name="extrapolate initial guess" type="bool" value="true"/>
      <Parameter name="nonlinear iteration initial guess extrapolation order" type="int" value="1"/>
      <Parameter name="restart tolerance relaxation factor" type="double" value="1.0"/>
      <Parameter name="restart tolerance relaxation factor damping" type="double" value="1.0"/>

      <Parameter name="timestep controller type" type="string" value="standard"/>
      <ParameterList name="timestep controller standard parameters">
        ...
      </ParameterList>

      <Parameter name="solver type" type="string" value="nka"/>
      <ParameterList name="nka parameters">
        ... 
      </ParameterList>
    </ParameterList>
  </ParameterList>




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

.. code-block:: xml

  <ParameterList name="BDF1"> <!-- parent list -->
    <Parameter name="timestep controller type" type="string" value="standard"/>
    <ParameterList name="timestep controller standard parameters">
      <Parameter name="min iterations" type="int" value="10"/>
      <Parameter name="max iterations" type="int" value="15"/>
      <Parameter name="time step increase factor" type="double" value="1.2"/>
      <Parameter name="time step reduction factor" type="double" value="0.5"/>
      <Parameter name="max time step" type="double" value="1e+9"/>
      <Parameter name="min time step" type="double" value="0.0"/>
    </ParameterList>
  </ParameterList>

In this example, the time step is increased by factor 1.2 when the nonlinear
solver converges in 10 or less iterations. 
The time step is not changed when the number of nonlinear iterations is
between 11 and 15.
The time step will be cut twice if the number of nonlinear iterations exceeds 15.




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
      happens when the current l2-error exceeds the initial l2-error.

 


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

.. code-block:: xml

  <Parameter name="solver type" type="string" value="JFNK"/>
  <ParameterList name="JFNK parameters">
    <Parameter name="typical solution value" type="double" value="1.0"/>

    <ParameterList name="JF matrix parameters">
      <Parameter name="finite difference epsilon" type="double" value="1.0e-8"/>
      <Parameter name="method for epsilon" type="string" value="Knoll-Keyes L2"/>
    </ParameterList>

    <ParameterList name="nonlinear solver">
      <Parameter name="nonlinear tolerance" type="double" value="1.0e-05"/>
      <Parameter name="diverged tolerance" type="double" value="1.0e+10"/>
      <Parameter name="limit iterations" type="int" value="20"/>
      <Parameter name="max divergent iterations" type="int" value="3"/>
    </ParameterList>

    <ParameterList name="linear operator">
      <Parameter name="iterative method" type="string" value="gmres"/>
      <ParameterList name="gmres parameters">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>





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

    * `"inner solver`" ``[solver-typed-spec]`` A Solver_, used at each step.




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

The interface to Trilinos NOX solver is as follows:

.. code-block:: xml

  <Parameter name="solver type" type="string" value="nox"/>
    <ParameterList name="nox parameters">
      <Parameter name="typical solution value" type="double" value="1.0"/>

      <ParameterList name="JF matrix parameters">
        <Parameter name="finite difference epsilon" type="double" value="1.0e-8"/>
        <Parameter name="method for epsilon" type="string" value="Knoll-Keyes L2"/>
      </ParameterList>

      <ParameterList name="nonlinear solver">
        <Parameter name="solver type" type="string" value="Newton"/>
        <ParameterList name="Newton parameters">
          ...
        </ParameterList>
      </ParameterList>

      <ParameterList name="linear operator">
        <Parameter name="iterative method" type="string" value="gmres"/>
        <ParameterList name="gmres parameters">
          ...
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>





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

.. code-block:: xml

  <ParameterList name="_PCG with HYPRE AMG">  <!-- parent list -->
  <ParameterList name="pcg parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual,make one iteration}"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>




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

.. code-block:: xml

  <ParameterList name="_GMRES with HYPRE AMG">  <!-- parent list -->
  <ParameterList name="gmres parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
    <Parameter name="size of Krylov space" type="int" value="10"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>
    <Parameter name="maximum size of deflation space" type="int" value="0"/>
    <Parameter name="preconditioning strategy`" type="string" value="left"/>
    <Parameter name="release Krylov vectors" type="bool" value="false"/>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>

  <!-- Alternative implementation
  <ParameterList name="belos gmres parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
    <Parameter name="size of Krylov space" type="int" value="10"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>
  </ParameterList-->
  </ParameterList>




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

.. code-block:: xml

  <ParameterList name="_NKA">  <!-- parent list -->
  <ParameterList name="nka parameters">
    <Parameter name="error tolerance" type="double" value="1e-12"/>
    <Parameter name="maximum number of iterations" type="int" value="400"/>
    <Parameter name="convergence criteria" type="Array(string)" value="{relative residual}"/>
    <Parameter name="overflow tolerance" type="double" value="3.0e+50"/>
    <Parameter name="max nka vectors" type="int" value="10"/>
    <Parameter name="nka vector tolerance" type="double" value="0.05"/>

    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>




Linear Solver: Amesos
---------------------
 Direct solvers via Trilinos.

Amesos library of Trilinos package provides interfaces to a few direct solvers.
List `"amesos parameters`" contains parameters that understood by this library.
These parameters may violate the camel-case convention employed by this spec.
Additional parameters are:

* `"solver name`" [string] declares name of one of the supported direct solvers. 
  Available options are `"klu`", `"superludist`", `"basker`", etc, see Amesos and 
  Amesos2 manuals for details. The default value is serial solver `"klu`".

* `"amesos version`" [int] specifies version of Amesos. Available options are 1 and 2.
  The default value is 1.

.. code-block:: xml

  <ParameterList name="_AMESOS KLU">  <!-- parent list -->
  <Parameter name="direct method" type="string" value="amesos"/>
  <ParameterList name="amesos parameters">
    <Parameter name="solver name" type="string" value="klu"/>
    <Parameter name="amesos version" type="int" value="1"/>
  </ParameterList>
  </ParameterList>




Linear Solver: Amesos
---------------------
 Direct solvers via Trilinos.

.. warning:: undocumented




Linear Solver: Belos (GMRES)
----------------------------
 Trilinos/Belos implementations of iterative methods.

Includes GMRES and other Belos methods.

.. warning:: undocumented





Preconditioners
===============
.. _Preconditioner:

 A base class for assembled preconditioners.

This sublist contains entries for various
preconditioners required by a simulation. At the moment, we support Trilinos multilevel
preconditioner, Hypre BoomerAMG preconditioner, ILU preconditioner, Hypre's Euclid ILU
preconditioner, and identity preconditioner. 

* `"preconditioning method`" [string] defines preconditioner algorithm.

* `"xxx parameters`" [list] provides parameters for the preconditioner specified 
  by parameter `"preconditioning method`".
 
.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="preconditioners">
    <ParameterList name="_TRILINOS ML">
      <Parameter name="preconditioning method" type="string" value="ml"/>
      <ParameterList name="ml parameters">
        ... 
      </ParameterList>
    </ParameterList>

    <ParameterList name="_HYPRE AMG">
      <Parameter name="preconditioning method" type="string" value="boomer amg"/>
      <ParameterList name="boomer amg parameters">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="_BLOCK ILU">
      <Parameter name="preconditioning method" type="string" value="block ilu"/>
      <ParameterList name="block ilu parameters">
        ...
      </ParameterList>
    </ParameterList>

    <ParameterList name="_DIAGONAL">
      <Parameter name="preconditioning method" type="string" value="diagonal"/>
    </ParameterList>
  </ParameterList>




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

The full list of relevant parameters is somewhat method-dependent, and is
documented extensively in the `Trilinos Documentation
<https://docs.trilinos.org/dev/packages/ifpack/doc/html/index.html>`_.

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
The algorithm was presented at SC99 and published in expanded 
form in the SIAM Journal on Scientific Computing. 
Scalability means that the factorization (setup) and application (triangular solve) timings remain
nearly constant when the global problem size is scaled in proportion to the number of processors.
As with all ILU preconditioning methods, the number of iterations is expected to increase with
global problem size.

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

   * `"cycles start period stop`" ``[Array(int)]`` **optional** The first entry
     is the start cycle, the second is the cycle period, and the third is the
     stop cycle or -1, in which case there is no stop cycle. A visualization
     dump is written at such cycles that satisfy cycle = start + n*period, for
     n=0,1,2,... and cycle < stop if stop != -1.0.

   * `"cycles start period stop 0`" ``[Array(int)]`` **optional** If multiple
     cycles start period stop parameters are needed, then use these parameters.
     If one with 0 is found, then one with 1 is looked for, etc, until the Nth
     one is not found.

   * `"cycles`" ``[Array(int)]`` **optional** An array of discrete cycles that
     at which a visualization dump is written.

   * `"times start period stop`" ``[Array(double)]`` **optional** The first
     entry is the start time, the second is the time period, and the third is
     the stop time or -1, in which case there is no stop time. A visualization
     dump is written at such times that satisfy time = start + n*period, for
     n=0,1,2,... and time < stop if stop != -1.0.

   * `"times start period stop units`" ``[string]`` **s** Units corresponding
     to this spec.  One of `"s`", `"min`", `"h`", `"d`", `"yr`", or `"noleap`"

   * `"times start period stop 0`" ``[Array(double)]`` **optional** If multiple
     start period stop parameters are needed, then use this these parameters
     with N=0,1,2.  If one with 0 is found, then one with 1 is looked for, etc,
     until the Nth one is not found.

   * `"times start period stop 0 units`" ``[string]`` **s** Units corresponding
     to this spec.  One of `"s`", `"min`", `"h`", `"d`", `"yr`", or `"noleap`"
     See above for continued integer listings.

   * `"times`" ``[Array(double)]`` **optional** An array of discrete times that
     at which a visualization dump shall be written.

   * `"times units`" ``[string]`` **s** Units corresponding to this spec.  One
     of `"s`", `"d`", `"yr`", or `"yr 365`"




Verbose Object
==============


This allows control of log-file verbosity for a wide variety of objects
and physics.

.. _verbose-object-spec:
.. admonition:: verbose-object-spec

   * `"verbosity level`" ``[string]`` **GLOBAL_VERBOSITY**, `"low`",
     `"medium`", `"high`", `"extreme`" The default is set by the global
     verbosity spec, (fix me!)  Typically, `"low`" prints out minimal
     information, `"medium`" prints out errors and overall high level
     information, `"high`" prints out basic debugging, and `"extreme`" prints
     out local debugging information.

   * `"write on rank`" ``[int]`` **0** VerboseObjects only write on a single
     rank -- by deafult the 0th rank.  However, sometimes it is useful for
     debugging to write from another rank due to a need for cleaner output or
     writing a specific cell/entity information.

   * `"output filename`" ``[string]`` **optional** Redirect this output to a
     specific file rather than writing to screen.  Note this will be done
     by-the-instance, so this may not catch as much as one might think.

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
 A base class for all functions of space and time.

Analytic, algabraic functions of space and time are used for a variety of
purposes, including boundary conditions, initial conditions, and independent
variables.

For initial conditions, functions are prescribed of space only, i.e.

.. math::
   u = f(x,y,z)

For boundary conditions and independent variables, functions are also a
function of time:

.. math::
   u = f(t,x,y,z)




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

.. _function-smooth-step-spec:
.. admonition:: function-smooth-step-spec

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

.. _function-polynomial-spec:
.. admonition:: function-polynomial-spec

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

.. _function-linear-spec:
.. admonition:: function-linear-spec

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
    <Parameter name="gradient" type="Array(double)" value="{1.0, 2.0, 3.0}"/>
    <Parameter name="x0" type="Array(double)" value="{2.0, 3.0, 1.0}"/>
  </ParameterList>


  

Separable Function
------------------
 FunctionSeparable: f(x,y) = f1(x)*f2(y)

A separable function is defined as the product of other functions such as

.. math::
  f(x_0, x_1,...,x_{{n-1}}) = f_1(x_0)\, f_2(x_1,...,x_{{n-1}})

where :math:`f_1` is defined by the `"function1`" sublist, and
:math:`f_2` by the `"function2`" sublist.

.. _function-separable-spec:
.. admonition:: function-separable-spec

   * `"function1`" ``[function-spec]`` :math:`f_1` in :math:`f(x) = f_1(x0) * f_2(x1...)`
   * `"function2`" ``[function-spec]`` :math:`f_2` in :math:`f(x) = f_1(x0) * f_2(x1...)`


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

.. _function-additive-spec:
.. admonition:: function-additive-spec

   * `"function1`" ``[function-spec]`` :math:`f_1` in :math:`f(x) = f_1(x) + f_2(x)`
   * `"function2`" ``[function-spec]`` :math:`f_2` in :math:`f(x) = f_1(x) + f_2(x)`

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

.. _function-multiplicative-spec:
.. admonition:: function-multiplicative-spec

   * `"function1`" ``[function-spec]`` :math:`f_1` in :math:`f(x) = f_1(x) * f_2(x)`
   * `"function2`" ``[function-spec]`` :math:`f_2` in :math:`f(x) = f_1(x) * f_2(x)`

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

.. _function-composition-spec:
.. admonition:: function-composition-spec

   * `"function1`" ``[function-spec]`` :math:`f_1` in :math:`f(x) = f_1(f_2(x))`
   * `"function2`" ``[function-spec]`` :math:`f_2` in :math:`f(x) = f_1(f_2(x))`


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

Given a two-dimensional array :math:`u_{i,j}`, :math:`f` is then defined by
bilinear interpolation on :math:`u_{i(x),j(y)}, u_{i(x)+1,j(y)}, u_{i(x),j(y)+1}, u_{i(x)+1,j(y)+1}`, 
if :math:`(x,y)` is in :math:`[x_0,x_n] \times [y_0,y_m]`, linear interpolation if one of :math:`x,y`
are out of those bounds, and constant at the corner value if both are out of
bounds.

.. _function-bilinear-spec:
.. admonition:: function-bilinear-spec

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

.. _function-distance-spec:
.. admonition:: function-distance-spec

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

.. _function-monomial-spec:
.. admonition:: function-monomial-spec

   * `"c`" ``[double]`` c in :math:`f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}`
   * `"x0`" ``[Array(double)]`` x0 in :math:`f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}`
   * `"exponents`" ``[Array(int)]`` p in :math:`f = c \prod_{j=0}^{n} (x_j - x_{0,j})^{p_j}`

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
 Provides access to many common mathematical functions.

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

.. _function-standard-math-spec:
.. admonition:: function-standard-math-spec

   * `"operator`" ``[string]`` specifies the name of a standard mathematical
     function.  Available options are:

     - trigonometric operators: `"cos`", `"sin`", `"tan`", `"acos`", `"asin`",
       `"atan`"
     - hyperbolic trig operators: `"cosh`", `"sinh`", `"tanh`"
     - power/log operators: `"pow`", `"exp`", `"log`", `"log10`", `"sqrt`",
     - integral operators: `"ceil`", `"floor`", `"mod`", 
     - `"abs`", `"fabs`", `"positive`" (0 for negative values), `"negative`" (0
       for positive values), `"heaviside`", `"sign`"

   * `"amplitude`" ``[double]`` specifies a multiplication factor `a` in
     formula `a f(x)`.  The multiplication factor is ignored by function
     `mod`. Default value is 1.

   * `"parameter`" ``[double]`` **1.0** specifies additional parameter `p` for
     math functions with two arguments. These functions are `"a pow(x[0], p)`"
     and `"a mod(x[0], p)`".  Alternative, scales the argument before
     application, for use in changing the period of trig functions.

   * `"shift`" ``[double]`` specifies a shift of the function argument. Default
     is 0.

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

Operators are discrete forms of linearized PDEs operators.
They form a layer between physical process kernels and solvers
and include accumulation, diffusion, advection, elasticity, reaction, 
and source operators.
The residual associated with an operator :math:`L_h` helps to 
understand the employed sign convention:

.. math::
  r = f - L_h u.

Operator represents a map from linear space X to linear space Y.  Typically,
this map is a linear map, and encapsulates much of the discretization involved
in moving from continuous to discrete equations. The spaces X and Y are described
by CompositeVectors (CV). A few maps X->Y are supported.

An operator provides an interface for applying both the forward and inverse
linear map (assuming the map is invertible).

Typically the Operator is never seen by the user; instead the user provides
input information for helper classes based on the continuous mathematical
operator and the desired discretization.  These helpers build the needed
``Operator``, which may include information from multiple helpers (i.e. in the
case of Jacobian Operators for a PDE).

However, one option may be provided by the user, which is related to dealing
with nearly singular operators:

* `"diagonal shift`" ``[double]`` **0.0** Adds a scalar shift to the diagonal
  of the ``Operator``, which can be useful if the ``Operator`` is singular or
  near-singular.

A PK decides how to bundle operators in a collection of operators.
For example, an advection-diffusion problem may benefit from using
a single operator that combines two operators representing diffusion and advection process.
Collection of operators must be used for implicit solvers and for building preconditioners.
In such a case, the collections acts as a single operator.

Operators use a few tools that are generic in nature and can be used independently by PKs. 
The list includes reconstruction and limiting algorithms. 


Schema
------

The operators use notion of schema to describe operator's abstract structure.
Old operators use a simple schema which is simply the list of geometric objects where
scalar degrees of freedom are defined.
New operators use a list to define location, type, and number of degrees of freedom.
In addition, the base of local stencil is either *face* or *cell*.
A rectangular operator needs two schemas do describe its domain (called `"schema domain`") 
and its range (called `"schema range`").
A square operator may use either two identical schema lists or a single list called `"schema`".

.. code-block:: xml

  <ParameterList name="pks operator name">  <!-- parent list-->
  <ParameterList name="schema domain">
    <Parameter name="base" type="string" value="cell"/>
    <Parameter name="location" type="Array(string)" value="{node, face}"/>
    <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
    <Parameter name="number" type="Array(int)" value="{2, 1}"/>
  </ParameterList>
  <ParameterList name="schema domain">
    <Parameter name="base" type="string" value="cell"/>
    <Parameter name="location" type="Array(string)" value="{node, face}"/>
    <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
    <Parameter name="number" type="Array(int)" value="{2, 1}"/>
  </ParameterList>
  </ParameterList>

This example describes a square operator with two degrees of freedom per mesh node and one
degree of freedom per mesh face. 
The face-based degree of freedom is the normal component of a vector field. 
Such set of degrees of freedom is used in the Bernardi-Raugel element for discretizing 
Stokes equations.
Parameter `"base`" indicates that local matrices are associated with mesh cells. 




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
  \nabla \cdot K k (\nabla (u + b g z))

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

Diffusion is the most frequently used operator. It employs the old schema.

* `"pks operator name`" [list] a PK specific name for the diffusion operator.

  * `"discretization primary`" [string] specifies an advanced discretization method that
    has useful properties under some a priori conditions on the mesh and/or permeability tensor.
    The available options are `"mfd: optimized for sparsity`", `"mfd: optimized for monotonicity`",
    `"mfd: default`", `"mfd: support operator`", `"mfd: two-point flux approximation`",
    `"fv: default`", and `"nlfv: default`".
    The first option is recommended for general meshes.
    The second option is recommended for orthogonal meshes and diagonal absolute 
    permeability tensor. 

  * `"discretization secondary`" [string] specifies the most robust discretization method
    that is used when the primary selection fails to satisfy all a priori conditions.
    Default value is equal to that for the primary discretization.

  * `"diffusion tensor`" [string] specifies additional properties of the diffusion tensor.
    It allows us to solve problems with non-symmetric but positive definite tensors. 
    Available options are *symmetric* (default) and *nonsymmetric*.

  * `"nonlinear coefficient`" [string] specifies a method for treating nonlinear diffusion
    coefficient, if any. Available options are `"none`", `"upwind: face`", `"divk: cell-face`" (default),
    `"divk: face`", `"standard: cell`", and `"divk: cell-face-twin`".
    Symmetry preserving methods are the divk-family of methods and the classical cell-centered
    method (`"standard: cell`"). The first part of the name indicates the base scheme.
    The second part (after the semi-column) indicates required components of the composite vector
    that must be provided by a physical PK.
    Default is `"none`".

  * `"schema`" [Array(string)] defines the operator stencil. It is a collection of 
    geometric objects. It equals to `"{cell}`" for finite volume schemes. 
    It is typically `"{face, cell}`" for mimetic discretizations.

  * `"preconditioner schema`" [Array(string)] defines the preconditioner stencil.
    It is needed only when the default assembling procedure is not desirable. 
    If skipped, the `"schema`" is used instead. 

  * `"gravity`" [bool] specifies if flow is driven also by the gravity.

  * `"gravity term discretization`" [string] selects a model for discretizing the 
    gravity term. Available options are `"hydraulic head`" [default] and `"finite volume`". 
    The first option starts with equation for the shifted solution, i.e. the hydraulic head,
    and derives gravity discretization by the reserve shifting.
    The second option is based on the divergence formula.

  * `"gravity magnitude`" [double] defined magnitude of the gravity vector.

  * `"Newton correction`" [string] specifies a model for correction (non-physical) terms 
    that must be added to the preconditioner. These terms approximate some Jacobian terms.
    Available options are `"true Jacobian`" and `"approximate Jacobian`".
    The FV scheme accepts only the first options. The othre schemes accept only the second option.

  * `"scaled constraint equation`" [bool] rescales flux continuity equations on mesh faces.
    These equations are divided by the nonlinear coefficient. This option allows us to 
    treat the case of zero nonlinear coefficient. At moment this feature does not work 
    with non-zero gravity term. Default is *false*.

  * `"constraint equation scaling cutoff"`" [double] specifies the cutoff value for
    applying rescaling strategy described above.  

  * `"consistent faces`" [list] may contain a `"preconditioner`" and
    `"linear operator`" list (see sections Preconditioners_ and LinearSolvers_
    respectively).  If these lists are provided, and the `"discretization
    primary`" is of type `"mfd: *`", then the diffusion method
    UpdateConsistentFaces() can be used.  This method, given a set of cell
    values, determines the faces constraints that satisfy the constraint
    equation in MFD by assembling and inverting the face-only system.  This is
    not currently used by any Amanzi PKs.

  * `"fracture`" [Array(string)] provides list of regions that defines a fracture network.
    This parameter is used only by the coupled flow PK.


Example:

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
    <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
    <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
    <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
    <Parameter name="gravity" type="bool" value="true"/>
    <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
    <Parameter name="gravity magnitude" type="double" value="9.81"/>
    <Parameter name="nonlinear coefficient" type="string" value="upwind: face"/>
    <Parameter name="Newton correction" type="string" value="true Jacobian"/>

    <ParameterList name="consistent faces">
      <ParameterList name="linear solver">
        ...
      </ParameterList>
      <ParameterList name="preconditioner">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>

This example creates a p-lambda system, i.e. the pressure is
discretized in mesh cells and on mesh faces. 
The preconditioner is defined on faces only, i.e. cell-based unknowns
are eliminated explicitly and the preconditioner is applied to the
Schur complement.






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

* `"use manifold flux`"  ``[bool]`` **false** Computes the flux using algorithms
  and data structures for manifolds or fracture networks. 





Additional options for MFD with the gravity term include:

* `"gravity term discretization`" ``[string]`` selects a model for discretizing the
  gravity term. Available options are `"hydraulic head`" (default) and `"finite volume`".
  The first option starts with equation for the shifted solution, i.e. the hydraulic
  head, and derives gravity discretization by the reserve shifting.
  The second option is based on the divergence formula.








PDE_Advection
-------------



A high-order advection operator may have different domain and range and therefore requires two schemas.
The structure of the new schema is described in the previous section.
A high-order advection operator has two terms in a weak formulation, corresponding to 
volume and surface integrals. These two terms are discretixed using two operators with
matrix of types *advection* and *flux*, respectively.


* `"pks operator name`" [list] a PK specific name for the advection operator.

  * `"method`" [string] defines a discretization method. The available option is `"dg modal`".

  * `"method order`" [int] defines method order. For example, the classical low-order finite 
    volume scheme is equivalent to DG of order 0.

  * `"matrix type`" [string] defines matrix type. The supported options are `"advection`"
    and `"flux`".

  * `"dg basis`" [string] defines bases for DG schemes. The available options are 
    `"regularized`" (recommended), `"normalized`", `"orthonormalized`", and `"natural`" 
    (not recommended).

  * `"gradient operator on test function`" [bool] defines place of the gradient operator.
    For integration by parts schemes, the gradient is transfered to a test function.
    This option is needed for discretizing volumetric integrals.

  * `"jump operator on test function`" [bool] defines place of the jump operator.
    For integration by parts schemes, the jump operator is applied to a test function.
    This option is needed for discretizing surface fluxes.

  * `"flux formula`" [string] defines type of the flux. The available options 
    are `"Rusanov`" (default), `"upwind`", `"downwind`", and `"NavierStokes`".

  * `"schema domain`" [list] defines a discretization schema for the operator domain.

  * `"schema range`" [list] defines a discretization schema for the operator range. 

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="method" type="string" value="dg modal"/>
    <Parameter name="method order" type="int" value="2"/>
    <Parameter name="flux formula" type="string" value="Rusanov"/>
    <Parameter name="matrix type" type="string" value="flux"/>
    <Parameter name="jump operator on test function" type="bool" value="true"/>

    <ParameterList name="schema domain">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{node, face}"/>
      <Parameter name="type" type="Array(string)" value="{scalar, normal component}"/>
      <Parameter name="number" type="Array(int)" value="{2, 1}"/>
    </ParameterList>
    <ParameterList name="schema range">
      <Parameter name="base" type="string" value="cell"/>
      <Parameter name="location" type="Array(string)" value="{cell}"/>
      <Parameter name="type" type="Array(string)" value="{scalar}"/>
      <Parameter name="number" type="Array(int)" value="{1}"/>
    </ParameterList>
  </ParameterList>

In this example, we construct an operator for volumetric integrals in a weak formulation
of advection problem.

The only low-order advection operator in Amanzi is the upwind operator. 
It employes the old schema.

.. code-block:: xml

  <ParameterList name="pks operator name">
    <Parameter name="base" type="string" value="face"/>
    <Parameter name="schema" type="Array(string)" value="{cell}"/>
    <Parameter name="method order" type="int" value="0"/>
    <Parameter name="matrix type" type="string" value="advection"/>
  </ParameterList>





``PDE_AdvectionUpwind`` assembles the discrete form of:

.. math::
  \nabla \cdot (q C)

which advects quantity :math:`C` with fluxes :math:`q`.

This is a simple, first-order donor-upwind scheme, and is recommended
for use in diffusion-dominated advection-diffusion equations.




Field Initializers
==================

Fields, also known by their underlying datatype, the CompositeVector,
can be initialized in a variety of ways.  These are used in a variety
of places as generic capability.

Function Initialization
-----------------------
.. _CompositeVectorFunction:

 Mesh Functions, evaluate a function on a mesh and stick the result in a vector.

CompositeVectorFunctions are ways of evaluating a piecewise function on a
mesh and sticking the result into a CompositeVector.

This is used in a variety of ways -- Initial Conditions, Boundary
Conditions, and Independent Variable Evaluators.

Typically any containing object of this spec is a list of these specs.  The
list is indexed by Region, and the regions (logically) should partition the
domain (or boundary of the domain in the case of BCs).

Each entry in that list is a:

.. _composite-vector-function-spec:
.. admonition:: composite-vector-function-spec

   ONE OF

   * `"region`" ``[string]`` Region on which this function is evaluated.

   OR

   * `"regions`" ``[Array(string)]`` List of regions on which this function is evaluated.

   END

   ONE OF

   * `"component`" ``[string]`` Mesh component to evaluate this on.  This is
     one of "cell", "face", "node", "edge", or "boundary_face". The last two
     may require additional conditions, such as a proper mesh initialization.
     The mask "*" could be used in place of the component name.

   OR

   * `"components`" ``[Array(string)]`` Mesh components to evaluate this on.
     This is some collection of "cell", "face", "node", "edge", and/or
     "boundary_face". The last two may require additional conditions, such as a
     proper mesh initialization.  The array with the single entry "*" could be
     used to initialize all existing components.

   END

   * `"function`" ``[function-typedinline-spec]`` The spec to provide the actual algebraic function.

 


Column File Initialization
--------------------------


Interpolate a depth-based, 1D column of data onto a mesh.  Values are
prescribed only to cells.  Expected is an HDF5 file in the format:

Depth coordinates z:

  /z[:] = (z_0, z_1, ... , z_n)
     z_0 = 0.0
     z_n >= max depth of mesh
     z_i > z_(i-1)

Function values u:

  /f[:] = (f_0(z_0), f_1(z_1), ..., f_n(z_n))

.. _column-initialization-spec
.. admonition:: column-initialization-spec

   * `"file`" ``[string]`` HDF5 filename
   * `"z header`" ``[string]`` name of the z-coordinate data: `z` above.  Depth
     coordinates (positive downward from the surface), [m]
   * `"f header`" ``[string]`` name of the function data: `f` above.

   ONE OF

   * `"surface sideset`" ``[string]`` Region on the surface domain from which
     to start to determine columns.

   OR

   * `"surface sidesets`" ``[Array(string)]`` Regions on the surface domain
     from which to start to determine columns.

   END



Exodus File Initialization
--------------------------

