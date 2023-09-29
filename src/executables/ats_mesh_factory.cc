/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Simple wrapper that takes a ParameterList and generates all needed meshes.
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "MeshLogicalFactory.hh"
#include "MeshColumn.hh"
#include "MeshSurfaceCell.hh"
#include "GeometricModel.hh"

#include "ats_mesh_factory.hh"

namespace ATS {
namespace Mesh {

using namespace Amanzi;

//
// Create a mesh from an ExodusII file
//
// Collective on comm
Teuchos::RCP<AmanziMesh::Mesh>
createMeshFromFile(const std::string& mesh_name,
                   Teuchos::ParameterList& mesh_plist,
                   const Comm_ptr_type& comm,
                   const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                   State& S,
                   VerboseObject& vo)
{
  Teuchos::ParameterList& mesh_file_plist = mesh_plist.sublist("read mesh file parameters");
  auto mesh_factory_plist = Teuchos::rcp(new Teuchos::ParameterList("mesh factory"));

  // partitioner
  std::string partitioner = mesh_plist.get<std::string>("partitioner", "zoltan_rcb");
  mesh_factory_plist->sublist("unstructured").sublist("expert").set("partitioner", partitioner);

  // vo
  if (mesh_plist.isSublist("verbose object"))
    mesh_factory_plist->set("verbose object", mesh_plist.sublist("verbose object"));

  // file name
  std::string file;
  if (mesh_file_plist.isParameter("file")) {
    file = mesh_file_plist.get<std::string>("file");
  } else {
    Errors::Message msg("\"read mesh file\" list missing \"file\" parameter.");
    Exceptions::amanzi_throw(msg);
  }

  // create the MSTK factory and mesh
  AmanziMesh::MeshFactory factory(comm, gm, mesh_factory_plist);
  auto mesh = factory.create(file);

  if (mesh != Teuchos::null) {
    // potentially build columns
    if (mesh_plist.isParameter("build columns from set")) {
      std::string regionname = mesh_plist.get<std::string>("build columns from set");
      mesh->build_columns(regionname);
    } else if (mesh_plist.get("build columns", true)) {
      mesh->build_columns();
    }

    // verify
    checkVerifyMesh(mesh_plist, mesh);
  }

  // check for deformable
  bool deformable = mesh_plist.get<bool>("deformable mesh", false);
  S.RegisterMesh(mesh_name, mesh, deformable);
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  Registered mesh \"" << mesh_name << "\"." << std::endl;
  }
  return mesh;
}


//
// Generate a structured mesh.
//
// Collective on comm
Teuchos::RCP<AmanziMesh::Mesh>
createMeshGenerated(const std::string& mesh_name,
                    Teuchos::ParameterList& mesh_plist,
                    const Comm_ptr_type& comm,
                    const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                    State& S,
                    VerboseObject& vo)
{
  Teuchos::ParameterList& mesh_generated_plist = mesh_plist.sublist("generate mesh parameters");
  auto mesh_factory_plist = Teuchos::rcp(new Teuchos::ParameterList("mesh factory"));

  // partitioner
  std::string partitioner = mesh_plist.get<std::string>("partitioner", "zoltan_rcb");
  mesh_factory_plist->sublist("unstructured").sublist("expert").set("partitioner", partitioner);

  // vo
  if (mesh_plist.isSublist("verbose object"))
    mesh_factory_plist->set("verbose object", mesh_plist.sublist("verbose object"));

  // create mesh
  AmanziMesh::MeshFactory factory(comm, gm, mesh_factory_plist);
  auto mesh = factory.create(mesh_generated_plist);

  if (mesh != Teuchos::null) {
    // build columns
    if (mesh_plist.isParameter("build columns from set")) {
      std::string regionname = mesh_plist.get<std::string>("build columns from set");
      mesh->build_columns(regionname);
    } else if (mesh_plist.get("build columns", false)) {
      mesh->build_columns();
    }

    // verify
    checkVerifyMesh(mesh_plist, mesh);
  }

  // check for deformable
  bool deformable = mesh_plist.get<bool>("deformable mesh", false);
  S.RegisterMesh(mesh_name, mesh, deformable);
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  Registered mesh \"" << mesh_name << "\"." << std::endl;
  }
  return mesh;
}


//
// Generate a logical mesh.
//
// Not collective (logical meshes are currently serial!)
Teuchos::RCP<AmanziMesh::Mesh>
createMeshLogical(const std::string& mesh_name,
                  Teuchos::ParameterList& mesh_plist,
                  const Comm_ptr_type& comm,
                  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                  State& S,
                  VerboseObject& vo)

{
  Teuchos::ParameterList& mesh_logical_plist = mesh_plist.sublist("logical mesh parameters");

  // -- from logical mesh file, has its own factory
  AmanziMesh::MeshLogicalFactory fac(comm, gm);

  Teuchos::RCP<AmanziMesh::Mesh> mesh = Teuchos::null;
  if (mesh_logical_plist.isParameter("read from file")) {
    auto filename = mesh_logical_plist.get<std::string>("read from file");
    auto my_list_in_other_file = Teuchos::getParametersFromXmlFile(filename);
    mesh = fac.Create(*my_list_in_other_file);
  } else {
    mesh = fac.Create(mesh_logical_plist);
  }

  if (mesh != Teuchos::null) { checkVerifyMesh(mesh_plist, mesh); }

  bool deformable = mesh_plist.get<bool>("deformable mesh", false);
  S.RegisterMesh(mesh_name, mesh, deformable);
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  Registered mesh \"" << mesh_name << "\"." << std::endl;
  }
  return mesh;
}


//
// Get an aliased mesh.
//
// Not collective.
Teuchos::RCP<const AmanziMesh::Mesh>
createMeshAliased(const std::string& mesh_name,
                  Teuchos::ParameterList& mesh_plist,
                  State& S,
                  VerboseObject& vo)
{
  Teuchos::ParameterList& alias_plist = mesh_plist.sublist("aliased parameters");

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = Teuchos::null;
  std::string target;
  if (alias_plist.isParameter("target")) {
    target = alias_plist.get<std::string>("target");

    // if both are domain set, construct the target from alias's index and
    // target's prefix
    KeyTriple ds_target, ds_alias;
    if (Keys::splitDomainSet(target, ds_target) && Keys::splitDomainSet(mesh_name, ds_alias)) {
      if (std::get<1>(ds_target) == "*") {
        target = Keys::getDomainInSet(std::get<0>(ds_target), std::get<1>(ds_alias));
      }
    }

    if (S.HasMesh(target)) {
      mesh = S.GetMesh(target);
    } else {
      Errors::Message msg("Aliased mesh \"");
      msg << mesh_name << "\" target mesh \"" << target << "\" does not exist in State."
          << "  Ensure it appears in the \"mesh\" list prior to the aliased mesh.";
      Exceptions::amanzi_throw(msg);
    }

  } else {
    Errors::Message msg("Aliased mesh \"");
    msg << mesh_name << "\" is missing parameter \"target\"";
    Exceptions::amanzi_throw(msg);
  }

  if (mesh != Teuchos::null) S.AliasMesh(target, mesh_name);
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  Aliased mesh \"" << mesh_name << "\" to \"" << target << "\"." << std::endl;
  }
  return mesh;
}


//
// Create a surface mesh by extracting from a volume mesh.
//
// Collective on _volume_ mesh communicator.
Teuchos::RCP<AmanziMesh::Mesh>
createMeshSurface(const std::string& mesh_name,
                  Teuchos::ParameterList& mesh_plist,
                  const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                  State& S,
                  VerboseObject& vo)
{
  Teuchos::ParameterList& mesh_surface_plist = mesh_plist.sublist("surface parameters");
  auto mesh_factory_plist = Teuchos::rcp(new Teuchos::ParameterList("mesh factory"));

  // vo
  if (mesh_plist.isSublist("verbose object"))
    mesh_factory_plist->set("verbose object", mesh_plist.sublist("verbose object"));

  Teuchos::RCP<AmanziMesh::Mesh> mesh = Teuchos::null;

  // get the parent mesh
  auto parent_name = mesh_surface_plist.get<std::string>("parent domain", "domain");
  if (Keys::isDomainSet(parent_name)) {
    KeyTriple dset_parent;
    Keys::splitDomainSet(parent_name, dset_parent);
    KeyTriple dset;
    Keys::splitDomainSet(mesh_name, dset);
    parent_name = Keys::getDomainInSet(std::get<0>(dset_parent), std::get<1>(dset));
  }
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  extracted from \"" << parent_name << "\"." << std::endl;
  }

  if (S.HasMesh(parent_name)) {
    auto parent = S.GetMesh(parent_name);
    auto comm = parent->get_comm();

    // get the regions
    std::vector<std::string> regions;
    if (mesh_surface_plist.isParameter("surface sideset name")) {
      regions.push_back(mesh_surface_plist.get<std::string>("surface sideset name"));
    } else if (mesh_surface_plist.isParameter("surface sideset names")) {
      regions =
        mesh_surface_plist.get<Teuchos::Array<std::string>>("surface sideset names").toVector();
    } else if (mesh_surface_plist.isParameter("region")) {
      regions.push_back(mesh_surface_plist.get<std::string>("region"));
    } else if (mesh_surface_plist.isParameter("regions")) {
      regions = mesh_surface_plist.get<Teuchos::Array<std::string>>("regions").toVector();
    } else {
      Errors::Message msg;
      msg << "Mesh \"" << mesh_name << "\" of type \"surface\" is missing parameter \"regions\".";
      Exceptions::amanzi_throw(msg);
    }

    // create the MSTK factory
    AmanziMesh::MeshFactory factory(comm, gm, mesh_factory_plist);
    Teuchos::RCP<AmanziMesh::Mesh> surface3D_mesh = Teuchos::null;

    // construct a 3D submanifold mesh if needed and the flattened surface mesh
    if (parent->manifold_dimension() == 3) {
      surface3D_mesh = factory.create(parent, regions, AmanziMesh::FACE, false, true, false);
      mesh = factory.create(parent, regions, AmanziMesh::FACE, true, true, false);
    } else {
      mesh = factory.create(parent, regions, AmanziMesh::CELL, true, true, false);
    }

    bool deformable = mesh_plist.get<bool>("deformable mesh", false);

    // register with state
    if (mesh != Teuchos::null) {
      checkVerifyMesh(mesh_plist, mesh);
      S.RegisterMesh(mesh_name, mesh, deformable);
      if (vo.os_OK(Teuchos::VERB_HIGH)) {
        *vo.os() << "  Registered mesh \"" << mesh_name << "\"." << std::endl;
      }

      std::string mesh3d_name = mesh_name + "_3d";
      if (surface3D_mesh != Teuchos::null) {
        S.RegisterMesh(mesh3d_name, surface3D_mesh, deformable);
        if (vo.os_OK(Teuchos::VERB_HIGH)) {
          *vo.os() << "  Registered mesh \"" << mesh3d_name << "\"." << std::endl;
        }
      } else {
        auto target = mesh_surface_plist.get<std::string>("parent domain", "domain");
        S.AliasMesh(target, mesh3d_name);
        if (vo.os_OK(Teuchos::VERB_HIGH)) {
          *vo.os() << "  Aliased mesh \"" << mesh3d_name << "\" to \"" << target << "\"."
                   << std::endl;
        }
      }
    }
  }
  return mesh;
}


//
// Create a mesh by extracting from another mesh.
//
// Collective on the _parent_ mesh communicator.
Teuchos::RCP<AmanziMesh::Mesh>
createMeshExtracted(const std::string& mesh_name,
                    Teuchos::ParameterList& mesh_plist,
                    const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                    State& S,
                    VerboseObject& vo)
{
  Teuchos::ParameterList& mesh_extracted_plist = mesh_plist.sublist("extracted parameters");
  auto mesh_factory_plist = Teuchos::rcp(new Teuchos::ParameterList("mesh factory"));

  // vo
  if (mesh_plist.isSublist("verbose object"))
    mesh_factory_plist->set("verbose object", mesh_plist.sublist("verbose object"));

  Teuchos::RCP<AmanziMesh::Mesh> mesh = Teuchos::null;

  // get the parent mesh
  auto parent_name = mesh_extracted_plist.get<std::string>("parent domain", "domain");
  if (Keys::isDomainSet(parent_name)) {
    KeyTriple dset_parent;
    Keys::splitDomainSet(parent_name, dset_parent);
    KeyTriple dset;
    Keys::splitDomainSet(mesh_name, dset);
    parent_name = Keys::getDomainInSet(std::get<0>(dset_parent), std::get<1>(dset));
  }
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  extracted from \"" << parent_name << "\"." << std::endl;
  }

  if (S.HasMesh(parent_name)) {
    auto parent = S.GetMesh(parent_name);
    auto comm = parent->get_comm();

    // get the regions
    std::vector<std::string> regions;
    if (mesh_extracted_plist.isParameter("region")) {
      regions.push_back(mesh_extracted_plist.get<std::string>("region"));
    } else if (mesh_extracted_plist.isParameter("regions")) {
      regions = mesh_extracted_plist.get<Teuchos::Array<std::string>>("regions").toVector();
    } else {
      Errors::Message msg;
      msg << "Mesh \"" << mesh_name << "\" of type \"extracted\" is missing parameter \"regions\".";
      Exceptions::amanzi_throw(msg);
    }

    // create the MSTK factory
    AmanziMesh::MeshFactory factory(comm, gm, mesh_factory_plist);

    // construct the extracted mesh
    mesh = factory.create(parent, regions, AmanziMesh::CELL, false, true, false);
    bool deformable = mesh_plist.get<bool>("deformable mesh", false);

    // register with state
    if (mesh != Teuchos::null) {
      checkVerifyMesh(mesh_plist, mesh);
      S.RegisterMesh(mesh_name, mesh, deformable);
      if (vo.os_OK(Teuchos::VERB_HIGH)) {
        *vo.os() << "  Registered mesh \"" << mesh_name << "\"." << std::endl;
      }
    }
  }
  return mesh;
}


//
// Create a mesh by extracting a column of cells into a ColumnMesh
//
// Not collective -- Column meshes are serial.
Teuchos::RCP<AmanziMesh::Mesh>
createMeshColumn(const std::string& mesh_name,
                 Teuchos::ParameterList& mesh_plist,
                 const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                 State& S,
                 VerboseObject& vo)
{
  Teuchos::ParameterList& mesh_column_plist = mesh_plist.sublist("column parameters");

  AmanziMesh::Entity_ID lid = mesh_column_plist.get<AmanziMesh::Entity_ID>("entity LID");
  auto parent_name = mesh_column_plist.get<std::string>("parent domain", "domain");
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  Constructing MeshColumn of name " << mesh_name << " with parent " << parent_name
             << std::endl;
  }
  auto parent = S.GetMesh(parent_name);
  auto parent_list = Teuchos::rcp(new Teuchos::ParameterList(*parent->parameter_list()));
  auto mesh = AmanziMesh::createColumnMesh(parent, lid, parent_list);
  bool deformable = mesh_plist.get<bool>("deformable mesh", false);

  // build columns and verify
  if (mesh != Teuchos::null) {
    if (mesh_plist.isParameter("build columns from set")) {
      std::string regionname = mesh_plist.get<std::string>("build columns from set");
      mesh->build_columns(regionname);
    } else if (mesh_plist.get("build columns", false)) {
      mesh->build_columns();
    }
    checkVerifyMesh(mesh_plist, mesh);
  }
  S.RegisterMesh(mesh_name, mesh, deformable);
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  based on column LID: " << lid << std::endl
             << "  Registered mesh \"" << mesh_name << "\"." << std::endl;
  }
  return mesh;
}


//
// Create a mesh of a single cell -- a column's surface mesh.
//
// Not collective -- Column meshes are serial.
Teuchos::RCP<AmanziMesh::Mesh>
createMeshColumnSurface(const std::string& mesh_name,
                        Teuchos::ParameterList& mesh_plist,
                        const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                        State& S,
                        VerboseObject& vo)
{
  Teuchos::ParameterList& mesh_column_surf_plist = mesh_plist.sublist("column surface parameters");

  std::string parent_name = mesh_column_surf_plist.get<std::string>("parent domain");
  if (Keys::isDomainSet(parent_name)) {
    KeyTriple dset_parent;
    Keys::splitDomainSet(parent_name, dset_parent);
    KeyTriple dset;
    Keys::splitDomainSet(mesh_name, dset);
    parent_name = Keys::getDomainInSet(std::get<0>(dset_parent), std::get<1>(dset));
  }

  auto parent = S.GetMesh(parent_name);
  std::string surface_region = mesh_column_surf_plist.get<std::string>("surface region", "surface");
  if (vo.os_OK(Teuchos::VERB_HIGH))
    *vo.os() << "  Constructing MeshSurfaceCell of name " << mesh_name << " with parent "
             << parent_name << std::endl;
  auto mesh = Teuchos::rcp(new AmanziMesh::MeshSurfaceCell(parent, surface_region));
  bool deformable = mesh_plist.get<bool>("deformable mesh", false);

  checkVerifyMesh(mesh_plist, mesh);
  S.RegisterMesh(mesh_name, mesh, deformable);
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "  Registered mesh \"" << mesh_name << "\"." << std::endl;
  }
  return mesh;
}


//
// Create a collection of meshes indexed over a domain set.
//
// Collective over the mesh from which these are extracted.
//
// Domain set meshes relate to a bunch of other meshes, and these roles can be
//   filled by the same or different meshes in different cases.
//
// 1. The indexing parent.  The indexing parent determines how many meshes are
//   in the domain set.  The domain sets might be indexed by regions on this
//   mesh, or by entities on this mesh.
//
// 2. The extracting parent.  Most domain set subdomain meshes are formed via
//   extraction.
//
// 3. The accumulating parent.  This is a parent (which may or may not even
//   exist) upon which subdomains can be mapped back into.  This allows better
//   visualization.
//
// In the simplest case, e.g. extracting subdomains from a volume mesh, all
// three are the same (volume) mesh.
//
// Other valid cases are the corresponding surface meshes, where the indexing
// parent might be the volume mesh, the extracting parent of each subdomain
// surface is the subdomain extracted from the volume mesh, and the
// accumulating parent is the global surface mesh, which itself was extracted
// from the global volume mesh.
//
// Alternatively, subgrid models may extrude from or "hang off of" a mesh.
// These might be indexed from a mesh, but have no extracting or accumulating
// parent.
//
//
// An Indexed Domain Set is a set of meshes, one per entity in an indexing mesh.
void
createDomainSetIndexed(const std::string& mesh_name_pristine,
                       Teuchos::ParameterList& mesh_plist,
                       const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                       State& S,
                       VerboseObject& vo)
{
  // strip a :* from the end of the domain set name if needed
  std::string delim(1, Keys::dset_delimiter);
  std::string mesh_name;
  if (Keys::ends_with(mesh_name_pristine, delim + "*")) {
    mesh_name = mesh_name_pristine.substr(0, mesh_name_pristine.length() - 2);
  } else {
    mesh_name = mesh_name_pristine;
  }

  Teuchos::ParameterList& ds_list = mesh_plist.sublist("domain set indexed parameters");

  // get the indexing info
  auto regions = ds_list.get<Teuchos::Array<std::string>>("regions").toVector();
  auto entity_kind = AmanziMesh::entity_kind(ds_list.get<std::string>("entity kind"));
  std::string indexing_parent_name = ds_list.get<std::string>("indexing parent domain", "domain");

  if (S.HasMesh(indexing_parent_name)) {
    auto indexing_parent_mesh = S.GetMesh(indexing_parent_name);

    // is there a reference mesh for visualization?
    bool is_reference_mesh = ds_list.isParameter("referencing parent domain");
    std::string reference_parent_name;
    Teuchos::RCP<const AmanziMesh::Mesh> reference_mesh = Teuchos::null;
    if (is_reference_mesh) {
      reference_parent_name = ds_list.get<std::string>("referencing parent domain");
      if (S.HasMesh(reference_parent_name)) reference_mesh = S.GetMesh(reference_parent_name);
    } else {
      reference_parent_name = "NONE";
    }

    // for each subdomain, create a referencing map, a map from subdomain to reference mesh
    std::vector<std::string> subdomains;
    std::vector<int> lids;
    std::map<std::string, Teuchos::RCP<const std::vector<int>>> reference_maps;

    // if aliased, we deal with domain sets specially
    std::string alias_target;

    // create the subdomains, indexed over entities
    for (const auto& region : regions) {
      AmanziMesh::Entity_ID_List region_ents;
      indexing_parent_mesh->get_set_entities(
        region, entity_kind, AmanziMesh::Parallel_type::OWNED, &region_ents);
      const auto& map = indexing_parent_mesh->map(entity_kind, false);

      for (const AmanziMesh::Entity_ID& lid : region_ents) {
        // subdomain name
        AmanziMesh::Entity_ID gid = map.GID(lid);
        std::string subdomain = std::to_string(gid);
        subdomains.push_back(subdomain);
        std::string full_subdomain_name = Keys::getDomainInSet(mesh_name, subdomain);

        // set up the parameter list
        Teuchos::ParameterList subdomain_list;
        if (ds_list.isSublist(full_subdomain_name)) {
          subdomain_list = ds_list.sublist(full_subdomain_name);
        } else {
          subdomain_list = ds_list.sublist(Keys::getDomainInSet(mesh_name, "*"));
        }
        subdomain_list.setName(full_subdomain_name);

        auto subdomain_mesh_type = subdomain_list.get<std::string>("mesh type");
        auto& subdomain_param_list = subdomain_list.sublist(subdomain_mesh_type + " parameters");

        if (!subdomain_param_list.isParameter("entity GID"))
          subdomain_param_list.set("entity GID", gid);
        if (!subdomain_param_list.isParameter("entity LID"))
          subdomain_param_list.set("entity LID", lid);
        if (!subdomain_param_list.isParameter("entity kind"))
          subdomain_param_list.set("entity kind", AmanziMesh::entity_kind_string(entity_kind));
        if (!subdomain_param_list.isParameter("parent domain"))
          subdomain_param_list.set("parent domain", indexing_parent_name);

        // construct
        auto subdomain_mesh =
          createMesh(subdomain_list, indexing_parent_mesh->get_comm(), gm, S, vo);

        // create maps to the reference mesh
        if (is_reference_mesh) {
          // construct map into the reference mesh
          if (subdomain_mesh_type == "extracted" || subdomain_mesh_type == "column") {
            reference_maps[full_subdomain_name] = AmanziMesh::createMapToParent(*subdomain_mesh);
          } else if (subdomain_mesh_type == "surface" || subdomain_mesh_type == "column surface") {
            AMANZI_ASSERT(reference_mesh != Teuchos::null);
            reference_maps[full_subdomain_name] =
              AmanziMesh::createMapSurfaceToSurface(*subdomain_mesh, *reference_mesh);
          } else if (subdomain_mesh_type == "aliased") {
            // use the reference map from the target mesh, but first we have to determine the target mesh name
            alias_target = subdomain_param_list.get<std::string>("target");
            KeyTriple dset;
            bool is_ds = Keys::splitDomainSet(alias_target, dset);
            if (is_ds) alias_target = std::get<0>(dset);
          } else {
            Errors::Message msg;
            msg << "Mesh \"" << mesh_name
                << "\" domain set cannot create reference map to mesh of type \""
                << subdomain_mesh_type << "\".";
            Exceptions::amanzi_throw(msg);
          }
        }
      }
    }

    // construct and register the domain set
    Teuchos::RCP<AmanziMesh::DomainSet> ds = Teuchos::null;
    if (is_reference_mesh) {
      if (alias_target.empty()) {
        ds = Teuchos::rcp(new AmanziMesh::DomainSet(
          mesh_name, indexing_parent_mesh, subdomains, reference_mesh, reference_maps));
      } else {
        auto ref_domain_set = S.GetDomainSet(alias_target);
        AMANZI_ASSERT(ref_domain_set->get_referencing_parent() == reference_mesh);
        auto reference_maps = ref_domain_set->get_subdomain_maps();

        // these maps are all indexed by the target name, update to the aliased name.
        std::map<std::string, Teuchos::RCP<const std::vector<int>>> new_reference_maps;
        for (const auto& key_val : reference_maps) {
          KeyTriple old_ds;
          Keys::splitDomainSet(key_val.first, old_ds);
          new_reference_maps[Keys::getDomainInSet(mesh_name, std::get<1>(old_ds))] = key_val.second;
        }
        ds = Teuchos::rcp(new AmanziMesh::DomainSet(
          mesh_name, indexing_parent_mesh, subdomains, reference_mesh, new_reference_maps));
      }
    } else {
      ds = Teuchos::rcp(new AmanziMesh::DomainSet(mesh_name, indexing_parent_mesh, subdomains));
    }
    S.RegisterDomainSet(mesh_name, ds);
  }
}


// Region-based Domain Set is a set of meshes, one per region.
void
createDomainSetRegions(const std::string& mesh_name_pristine,
                       Teuchos::ParameterList& mesh_plist,
                       const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
                       State& S,
                       VerboseObject& vo)
{
  // strip a :* from the end of the domain set name if needed
  std::string delim(1, Keys::dset_delimiter);
  std::string mesh_name;
  if (Keys::ends_with(mesh_name_pristine, delim + "*")) {
    mesh_name = mesh_name_pristine.substr(0, mesh_name_pristine.length() - 2);
  } else {
    mesh_name = mesh_name_pristine;
  }

  Teuchos::ParameterList& ds_list = mesh_plist.sublist("domain set regions parameters");

  // get the indexing info
  auto regions = ds_list.get<Teuchos::Array<std::string>>("regions").toVector();
  std::string indexing_parent_name = ds_list.get<std::string>("indexing parent domain", "domain");

  if (S.HasMesh(indexing_parent_name)) {
    auto indexing_parent_mesh = S.GetMesh(indexing_parent_name);

    // is there a reference mesh for visualization?
    bool is_reference_mesh = ds_list.isParameter("referencing parent domain");
    std::string reference_parent_name;
    Teuchos::RCP<const AmanziMesh::Mesh> reference_mesh = Teuchos::null;
    if (is_reference_mesh) {
      reference_parent_name = ds_list.get<std::string>("referencing parent domain");
      if (S.HasMesh(reference_parent_name)) reference_mesh = S.GetMesh(reference_parent_name);
    } else {
      reference_parent_name = "NONE";
    }

    // for each subdomain, create a referencing map, a map from subdomain to reference mesh
    std::vector<std::string> subdomains;
    std::vector<int> lids;
    std::map<std::string, Teuchos::RCP<const std::vector<int>>> reference_maps;

    // create the subdomains, indexed over entities
    for (const auto& subdomain : regions) {
      std::string full_subdomain_name = Keys::getDomainInSet(mesh_name, subdomain);

      // set up the parameter list
      Teuchos::ParameterList subdomain_list;
      if (ds_list.isSublist(full_subdomain_name)) {
        subdomain_list = ds_list.sublist(full_subdomain_name);
      } else {
        subdomain_list = ds_list.sublist(Keys::getDomainInSet(mesh_name, "*"));
      }
      subdomain_list.setName(full_subdomain_name);

      auto subdomain_mesh_type = subdomain_list.get<std::string>("mesh type");
      auto& subdomain_param_list = subdomain_list.sublist(subdomain_mesh_type + " parameters");

      if (!subdomain_param_list.isParameter("parent domain"))
        subdomain_param_list.set("parent domain", indexing_parent_name);
      if (!subdomain_param_list.isParameter("region"))
        subdomain_param_list.set("region", subdomain);

      // construct
      auto subdomain_mesh = createMesh(subdomain_list, indexing_parent_mesh->get_comm(), gm, S, vo);

      if (subdomain_mesh != Teuchos::null) {
        subdomains.push_back(subdomain);

        // create maps to the reference mesh
        if (is_reference_mesh) {
          // construct map into the reference mesh
          if (subdomain_mesh_type == "extracted" || subdomain_mesh_type == "column") {
            reference_maps[full_subdomain_name] = createMapToParent(*subdomain_mesh);
          } else if (subdomain_mesh_type == "surface") {
            AMANZI_ASSERT(reference_mesh != Teuchos::null);
            reference_maps[full_subdomain_name] =
              createMapSurfaceToSurface(*subdomain_mesh, *reference_mesh);
          }
        }
      }
    }

    // construct and register the domain set
    Teuchos::RCP<AmanziMesh::DomainSet> ds = Teuchos::null;
    if (is_reference_mesh) {
      ds = Teuchos::rcp(new AmanziMesh::DomainSet(
        mesh_name, indexing_parent_mesh, subdomains, reference_mesh, reference_maps));
    } else {
      ds = Teuchos::rcp(new AmanziMesh::DomainSet(mesh_name, indexing_parent_mesh, subdomains));
    }
    S.RegisterDomainSet(mesh_name, ds);
  }
}

Teuchos::RCP<const Amanzi::AmanziMesh::Mesh>
createMesh(Teuchos::ParameterList& mesh_plist,
           const Amanzi::Comm_ptr_type& comm,
           const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
           Amanzi::State& S,
           Amanzi::VerboseObject& vo)
{
  auto mesh_type = mesh_plist.get<std::string>("mesh type");
  auto mesh_name = Keys::cleanPListName(mesh_plist.name());

  auto tab = vo.getOSTab();
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    *vo.os() << "Creating mesh \"" << mesh_name << "\" of type \"" << mesh_type << "\"."
             << std::endl;
  }

  if (mesh_type == "read mesh file") {
    return createMeshFromFile(mesh_name, mesh_plist, comm, gm, S, vo);
  } else if (mesh_type == "generate mesh") {
    return createMeshGenerated(mesh_name, mesh_plist, comm, gm, S, vo);
  } else if (mesh_type == "logical mesh") {
    return createMeshLogical(mesh_name, mesh_plist, comm, gm, S, vo);
  } else if (mesh_type == "aliased") {
    return createMeshAliased(mesh_name, mesh_plist, S, vo);
  } else if (mesh_type == "surface") {
    return createMeshSurface(mesh_name, mesh_plist, gm, S, vo);
  } else if (mesh_type == "extracted") {
    return createMeshExtracted(mesh_name, mesh_plist, gm, S, vo);
  } else if (mesh_type == "column") {
    return createMeshColumn(mesh_name, mesh_plist, gm, S, vo);
  } else if (mesh_type == "column surface") {
    return createMeshColumnSurface(mesh_name, mesh_plist, gm, S, vo);
  } else if (mesh_type == "domain set indexed") {
    createDomainSetIndexed(mesh_name, mesh_plist, gm, S, vo);
  } else if (mesh_type == "domain set regions") {
    createDomainSetRegions(mesh_name, mesh_plist, gm, S, vo);
  } else {
    Errors::Message msg;
    msg << "ATS Mesh Factory: unknown \"mesh type\" parameter \"" << mesh_type << "\" in mesh \""
        << mesh_name << "\".";
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}


bool
checkVerifyMesh(Teuchos::ParameterList& mesh_plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  // mesh verification
  AMANZI_ASSERT(!mesh.is_null());
  bool verify = mesh_plist.get<bool>("verify mesh", false);
  if (verify) {
    int num_procs = mesh->get_comm()->NumProc();
    int rank = mesh->get_comm()->MyPID();

    if (rank == 0) std::cout << "Verifying mesh with Mesh Audit..." << std::endl;
    if (num_procs == 1) {
      MeshAudit mesh_auditor(mesh);
      int status = mesh_auditor.Verify();
      if (status == 0) {
        std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
      } else {
        Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
        Exceptions::amanzi_throw(msg);
        return false;
      }
    } else {
      std::ostringstream ofile;
      ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
      std::ofstream ofs(ofile.str().c_str());
      if (rank == 0)
        std::cout << "Writing Mesh Audit output to " << ofile.str() << ", etc." << std::endl;

      int ierr = 0, aerr = 0;
      MeshAudit mesh_auditor(mesh, ofs);
      int status = mesh_auditor.Verify(); // check the mesh
      if (status != 0) ierr = 1;

      mesh->get_comm()->SumAll(&ierr, &aerr, 1);
      if (aerr == 0) {
        if (mesh->get_comm()->MyPID() == 0)
          std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
      } else {
        Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
        Exceptions::amanzi_throw(msg);
        return false;
      }
    }
  } // if verify
  return true;
}


void
createMeshes(Teuchos::ParameterList& global_list,
             const Comm_ptr_type& comm,
             const Teuchos::RCP<AmanziGeometry::GeometricModel>& gm,
             State& S)
{
  Teuchos::ParameterList& meshes_list = global_list.sublist("mesh");
  VerboseObject vo(comm, "ATS Mesh Factory", meshes_list);

  // always try to do the domain mesh first
  if (meshes_list.isSublist("domain")) {
    createMesh(meshes_list.sublist("domain"), comm, gm, S, vo);
  }

  // always try to do the surface mesh second
  if (meshes_list.isSublist("surface")) {
    createMesh(meshes_list.sublist("surface"), comm, gm, S, vo);
  }

  // now do the rest
  for (auto sublist : meshes_list) {
    if (sublist.first != "domain" && sublist.first != "surface" &&
        sublist.first != "verbose object" && meshes_list.isSublist(sublist.first)) {
      createMesh(meshes_list.sublist(sublist.first), comm, gm, S, vo);
    }
  }
}

} // namespace Mesh
} // namespace ATS
