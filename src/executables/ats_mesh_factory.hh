/*
  Copyright 2010-202x held jointly by participating institutions.
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
     - `"extracted`" See `Extracted Mesh`_.
     - `"domain set indexed`" See `Domain Set Meshes`_.
     - `"domain set regions`" See `Domain Set Meshes`_.
     - `"column`" See `Column Meshes`_.
     - `"column surface`" See `Column Surface Meshes`_.

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

`"mesh type`" = `"generate mesh`".

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

`"mesh type`" = `"read mesh file`".

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

`"mesh type`" = `"logical`".

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

`"mesh type`" = `"surface`".

.. _mesh-surface-spec:
.. admonition:: mesh-surface-spec

   * `"parent domain`" ``[string]`` **domain** Parent mesh's name.

   ONE OF

   * `"surface sideset name`" ``[string]`` The :doc:`region` name containing all surface faces.

   OR

   * `"surface sideset names`" ``[Array(string)]`` A list of :doc:`region` names containing the surface faces.

   END

   * `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
   * `"export mesh to file`" ``[string]`` **optional** Export the lifted
     surface mesh to this filename.
   * `"create subcommunicator`" ``[bool]`` **false** If false, the communicator
     of this mesh is the same as the parent mesh.  If true, the communicator of
     this mesh is the subset of the parent mesh comm that has entries on the
     surface.

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

Extracted Mesh
==============

A mesh is created by lifting a subset of entities from a parent mesh.  Locality
is preserved, so all local entities in this mesh have parents whose entities
are local on the parent mesh, so that no communication is ever done when
passing info between an parent mesh and an extracted mesh.

`"mesh type`" = `"extracted`".

.. _mesh-extracted-spec:
.. admonition:: mesh-extracted-spec

   * `"parent domain`" ``[string]`` **domain** Parent mesh's name.

   ONE OF

   * `"region`" ``[string]`` The Region_ name containing all surface faces.

   OR

   * `"regions`" ``[Array(string)]`` A list of Region_ names containing the surface faces.

   END

   * `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
   * `"export mesh to file`" ``[string]`` **optional** Export the lifted
     surface mesh to this filename.
   * `"create subcommunicator`" ``[bool]`` **false** If false, the communicator
     of this mesh is the same as the parent mesh.  If true, the communicator of
     this mesh is the subset of the parent mesh comm that has entries on the
     surface.
     

Aliased Mesh
============

Aliased domains are simply domains that share a mesh with another domain.  For
instance, one might find it useful to define both a "surface water" and a
"snow" domain that share a common "surface" mesh.  In that case, typically the
"surface" domain would point to the "surface" mesh, and the "snow" domain would
be an "aliased" domain whose target is the "surface" mesh.

`"mesh type`" = `"aliased`".

.. _mesh-aliased-spec:
.. admonition:: mesh-aliased-spec

   * `"target`" ``[string]`` Mesh that this alias points to.


Domain Set Meshes
=================

A collection of meshes formed by associating a new mesh with each entity of a
region or set of indices.  This includes generating a 1D column for each
surface face of a semi-structured subsurface mesh, or for hanging logical
meshes off of each surface cell as a subgrid model, etc.

The domain set meshes are then named `"MESH_NAME:X"` for each X, which can be a
local ID of an entity (in the case of `"domain set indexed`") or a region name
(in the case of `"domain set regions`").


`"mesh type`" = `"domain set indexed`".

.. _mesh-domain-set-indexed-spec:
.. admonition:: mesh-domain-set-indexed-spec

   * `"regions`" ``[Array(string)]`` Regions from which indices are created.
   * `"entity kind`" ``[string]`` One of `"cell`", `"face`", etc.  Entity of the
     region (usually `"cell`") on which each mesh will be associated.
   * `"indexing parent domain`" ``[string]`` **domain** Mesh which includes the above region.
   * `"referencing parent domain`" ``[string]`` **optional** Mesh from which
     the entities of the mesh will be extracted.  For instance, columns may be
     indexed from a surface mesh and referenced from the volume mesh below that
     surface.

Note, additionally, there must be a sublist of the name of the domain set,
which itself is a :doc:`mesh-typed-spec <mesh>`, but may be missing some info
(e.g. `"entity LID`") that is filled in by this index.

Example:

.. code-block:: xml

    <ParameterList name="column:*" type="ParameterList">
      <Parameter name="mesh type" type="string" value="domain set indexed" />
      <ParameterList name="domain set indexed parameters" type="ParameterList">
        <Parameter name="indexing parent domain" type="string" value="surface" />
        <Parameter name="entity kind" type="string" value="cell" />
        <Parameter name="referencing parent domain" type="string" value="domain" />
        <Parameter name="regions" type="Array(string)" value="{surface}" />
        <ParameterList name="column:*" type="ParameterList">
          <Parameter name="mesh type" type="string" value="column" />
          <ParameterList name="column parameters" type="ParameterList">
            <Parameter name="parent domain" type="string" value="domain" />
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>


`"mesh type`" = `"domain set regions`".

.. _mesh-domain-set-regions-spec:
.. admonition:: mesh-domain-set-regions-spec

   * `"regions`" ``[Array(string)]`` Regions from which indices are created.
   * `"entity kind`" ``[string]`` One of `"cell`", `"face`", etc.  Entity of the
     region (usually `"cell`") on which each mesh will be associated.
   * `"indexing parent domain`" ``[string]`` **domain** Mesh which includes the above region.
   * `"referencing parent domain`" ``[string]`` **optional** Mesh from which
     the entities of the mesh will be extracted.  For instance, columns may be
     indexed from a surface mesh and referenced from the volume mesh below that
     surface.

Note, additionally, there must be a sublist of the name of the domain set,
which itself is a :doc:`mesh-typed-spec <mesh>`, but may be missing some info
(e.g. `"region`") that is filled in by this domain set.

Example: the below example shows how to extract two subdomains, making them
each a proper mesh whose communicators only live where they have cells, thereby
decomposing the domain mesh into two subdomains.

.. code-block:: xml

    <ParameterList name="watershed:*" type="ParameterList">
      <Parameter name="mesh type" type="string" value="domain set regions" />
      <ParameterList name="domain set regions parameters" type="ParameterList">
        <Parameter name="indexing parent domain" type="string" value="domain" />
        <Parameter name="entity kind" type="string" value="cell" />
        <Parameter name="referencing parent domain" type="string" value="domain" />
        <Parameter name="regions" type="Array(string)" value="{upstream, downstream}" />
        <ParameterList name="watershed:*" type="ParameterList">
          <Parameter name="mesh type" type="string" value="extracted" />
          <ParameterList name="extracted parameters" type="ParameterList">
            <Parameter name="parent domain" type="string" value="domain" />
            <Parameter name="create subcommunicator" type="bool" value="true" />
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>


Column Meshes
=============

.. warning::
   Note these are rarely if ever created manually by a user.  Instead use
   `Domain Set Meshes`_, which generate a column mesh spec for every face
   of a set.

`"mesh type`" = `"column`".

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


Column Surface Meshes
=====================

.. warning::
   Note these are rarely if ever created manually by a user.  Instead use
   `Domain Set Meshes`_, which generate a column surface mesh spec for every face
   of a set.

`"mesh type`" = `"column surface`".

.. _mesh-column-surface-spec:
.. admonition:: mesh-column-surface-spec

   * `"parent domain`" ``[string]`` The name of the 3D mesh from which columns are generated.
     Note that the `"build columns from set`" parameter must be set in that mesh.
   * `"surface region`" ``[string]`` Region of the surface of the parent mesh.
   * `"verify mesh`" ``[bool]`` **false** Verify validity of surface mesh.
   * `"deformable mesh`" ``[bool]`` **false**  Used for deformation PKs to allow non-const access.


*/

#ifndef ATS_MESH_FACTORY_HH_
#define ATS_MESH_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "State.hh"
#include "VerboseObject.hh"


namespace ATS {
namespace Mesh {


//
// Helper function
//
bool
checkVerifyMesh(Teuchos::ParameterList& mesh_plist,
                Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh);

//
// Create mesh for each type
//
Teuchos::RCP<Amanzi::AmanziMesh::Mesh>
createMeshFromFile(const std::string& mesh_name,
                   const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                   const Amanzi::Comm_ptr_type& comm,
                   const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                   Amanzi::State& S,
                   Amanzi::VerboseObject& vo);

Teuchos::RCP<Amanzi::AmanziMesh::Mesh>
createMeshGenerated(const std::string& mesh_name,
                    const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                    const Amanzi::Comm_ptr_type& comm,
                    const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                    Amanzi::State& S,
                    Amanzi::VerboseObject& vo);

Teuchos::RCP<Amanzi::AmanziMesh::Mesh>
createMeshLogical(const std::string& mesh_name,
                  const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                  const Amanzi::Comm_ptr_type& comm,
                  const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                  Amanzi::State& S,
                  Amanzi::VerboseObject& vo);

Teuchos::RCP<const Amanzi::AmanziMesh::Mesh>
createMeshAliased(const std::string& mesh_name,
                  const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                  Amanzi::State& S,
                  Amanzi::VerboseObject& vo);

Teuchos::RCP<Amanzi::AmanziMesh::Mesh>
createMeshSurface(const std::string& mesh_name,
                  const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                  const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                  Amanzi::State& S,
                  Amanzi::VerboseObject& vo);

Teuchos::RCP<Amanzi::AmanziMesh::Mesh>
createMeshExtracted(const std::string& mesh_name,
                    const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                    const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                    Amanzi::State& S,
                    Amanzi::VerboseObject& vo);


Teuchos::RCP<Amanzi::AmanziMesh::Mesh>
createMeshColumn(const std::string& mesh_name,
                 const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                 const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                 Amanzi::State& S,
                 Amanzi::VerboseObject& vo);

Teuchos::RCP<Amanzi::AmanziMesh::Mesh>
createMeshColumnSurface(const std::string& mesh_name,
                        const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                        const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                        Amanzi::State& S,
                        Amanzi::VerboseObject& vo);

void
createDomainSetIndexed(const std::string& mesh_name_pristine,
                       const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                       const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                       Amanzi::State& S,
                       Amanzi::VerboseObject& vo);

void
createDomainSetRegions(const std::string& mesh_name_pristine,
                       const Teuchos::RCP<Teuchos::ParameterList>& mesh_plist,
                       const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
                       Amanzi::State& S,
                       Amanzi::VerboseObject& vo);

Teuchos::RCP<const Amanzi::AmanziMesh::Mesh>
createMesh(const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Amanzi::Comm_ptr_type& comm,
           const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
           Amanzi::State& s,
           Amanzi::VerboseObject& vo);

void
createMeshes(const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Amanzi::Comm_ptr_type& comm,
             const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
             Amanzi::State& s);

void
setDefaultParameters(Teuchos::ParameterList& plist, const Amanzi::VerboseObject& vo);

} // namespace Mesh
} // namespace ATS

#endif
