/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A list of mesh objects and their domain names.

/*!

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
    * `"format`" ``[string]`` format of pre-generated mesh file. One of:
    
      - `"MSTK`"
      - `"Exodus II`"

Example:

.. code-block:: xml

   <ParameterList name="mesh">
     <ParameterList name="domain">
       <Parameter name="mesh type" type="string" value="read mesh file"/>
       <ParameterList name="read mesh file parameters">
         <Parameter name="file" type="string" value="mesh_filename.exo"/>
         <Parameter name="format" type="string" value="Exodus II"/>
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

*/

#ifndef ATS_MESH_FACTORY_HH_
#define ATS_MESH_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "State.hh"


namespace ATS {

bool
checkVerifyMesh(Teuchos::ParameterList& mesh_plist,
                Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh);

void
createMesh(Teuchos::ParameterList& plist,
           const Amanzi::Comm_ptr_type& comm,
           const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
           Amanzi::State& s);


void
createMeshes(Teuchos::ParameterList& plist,
             const Amanzi::Comm_ptr_type& comm,             
             const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
             Amanzi::State& s);


} // namespace ATS

#endif
