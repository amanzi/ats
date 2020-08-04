/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Simple wrapper that takes a ParameterList and generates all needed meshes.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

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

void
createMesh(Teuchos::ParameterList& mesh_plist,
           const Amanzi::Comm_ptr_type& comm,           
           const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
           Amanzi::State& S)
{
  auto mesh_type = mesh_plist.get<std::string>("mesh type");
  auto mesh_name = Amanzi::Keys::cleanPListName(mesh_plist.name());
                   
  if (mesh_type == "read mesh file") {
    //Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    //prefs.clear();
    //prefs.push_back(Amanzi::AmanziMesh::MSTK);
    
    // from file
    auto read_params = Teuchos::rcp(new Teuchos::ParameterList(mesh_plist.sublist("read mesh file parameters")));
    
    // file name
    std::string file;
    if (read_params -> isParameter("file")) {
      file = read_params -> get<std::string>("file");
    } else {
      Errors::Message msg("\"read mesh file\" list missing \"file\" parameter.");
      Exceptions::amanzi_throw(msg);
    }

    // file format
    std::string format;
    if (read_params -> isParameter("format")) {
      // Is the format one that we can read?
      format = read_params -> get<std::string>("format");
      if (format != "Exodus II") {
        Errors::Message msg;
        msg << "\"read mesh file\" parameter \"format\" with value \"" << format
            << "\" not understood: valid formats are: \"Exodus II\".";
        Exceptions::amanzi_throw(msg);
      }
    } else {
      Errors::Message msg("\"read mesh file\" parameter \"format\" missing.");
      Exceptions::amanzi_throw(msg);
    }

    // create the MSTK factory
    Amanzi::AmanziMesh::MeshFactory factory(comm, gm, read_params);
    
    auto mesh = factory.create(file);

    if (mesh_plist.isParameter("build columns from set")) {
      std::string regionname = mesh_plist.get<std::string>("build columns from set");
      mesh->build_columns(regionname);
    }

    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(mesh_name, mesh, deformable);
    
  } else if (mesh_type == "generate mesh") {
    // create the MSTK factory
    // Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    // prefs.clear();
    // prefs.push_back(Amanzi::AmanziMesh::MSTK);

    // generated mesh

    auto genmesh_params = Teuchos::rcp(new Teuchos::ParameterList(mesh_plist.sublist("generate mesh parameters")));

    Amanzi::AmanziMesh::MeshFactory factory(comm, gm, genmesh_params);
    
    auto mesh = factory.create(*genmesh_params);
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    if (mesh_plist.isParameter("build columns from set")) {
      std::string regionname = mesh_plist.get<std::string>("build columns from set");
      mesh->build_columns(regionname);
    }

    checkVerifyMesh(mesh_plist, mesh);
    std::string name = mesh_plist.name();
    S.RegisterMesh(mesh_name, mesh, deformable);

  } else if (mesh_type == "logical mesh") {
    // -- from logical mesh file
    Amanzi::AmanziMesh::MeshLogicalFactory fac(comm, gm);
    auto log_mesh_params = mesh_plist.sublist("logical mesh parameters");
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
    if (log_mesh_params.isParameter("read from file")) {
      auto filename = log_mesh_params.get<std::string>("read from file");
      auto my_list_in_other_file = Teuchos::getParametersFromXmlFile(filename);
      mesh = fac.Create(*my_list_in_other_file);
    } else {
      mesh = fac.Create(log_mesh_params);
    }
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(mesh_name, mesh, deformable);

  } else if (mesh_type == "aliased") {
    S.AliasMesh(mesh_plist.sublist("aliased parameters").get<std::string>("alias"), mesh_name);

  } else if (mesh_type == "surface") {
    Teuchos::ParameterList& surface_plist = mesh_plist.sublist("surface parameters");
    std::vector<std::string> setnames;
    if (surface_plist.isParameter("surface sideset name")) {
      setnames.push_back(surface_plist.get<std::string>("surface sideset name"));
    } else if (surface_plist.isParameter("surface sideset names")) {
      setnames = surface_plist.get<Teuchos::Array<std::string> >("surface sideset names").toVector();
    } else {
      Errors::Message message("Surface mesh sublist missing parameter \"surface sideset names\".");
      Exceptions::amanzi_throw(message);
    }

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface3D_mesh = Teuchos::null;
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface_mesh = Teuchos::null;

    auto surface_plist_rcp = Teuchos::rcp(new Teuchos::ParameterList(surface_plist));
    // create the MSTK factory
    Amanzi::AmanziMesh::MeshFactory factory(comm, gm, surface_plist_rcp);
    // Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    // prefs.clear();
    // prefs.push_back(Amanzi::AmanziMesh::MSTK);

    auto parent = S.GetMesh(surface_plist.get<std::string>("parent domain", "domain"));
    if (parent->manifold_dimension() == 3) {
      surface3D_mesh = factory.create(parent,setnames,Amanzi::AmanziMesh::FACE,false,true,false);
      surface_mesh = factory.create(parent,setnames,Amanzi::AmanziMesh::FACE,true,true,false);
    } else {
      surface_mesh = factory.create(parent,setnames,Amanzi::AmanziMesh::CELL,true,true,false);
    }
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    if (surface3D_mesh.get()) {
      S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name())+"_3d", surface3D_mesh, deformable);
    } else {
      S.AliasMesh(surface_plist.get<std::string>("parent domain", "domain"),
                  Amanzi::Keys::cleanPListName(mesh_plist.name())+"_3d");
    }
    checkVerifyMesh(mesh_plist, surface_mesh);
    S.RegisterMesh(mesh_name, surface_mesh, deformable);

  } else if (mesh_type == "column") {
    Teuchos::ParameterList& column_list = mesh_plist.sublist("column parameters");
    Amanzi::AmanziMesh::Entity_ID lid = column_list.get<Amanzi::AmanziMesh::Entity_ID>("entity LID");
    auto parent = S.GetMesh(column_list.get<std::string>("column domain", "domain"));
    auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(parent, lid));
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(mesh_name, mesh, deformable);

  } else if (mesh_type == "column surface") {
    Teuchos::ParameterList& column_list = mesh_plist.sublist("column surface parameters");
    std::string surface_setname = column_list.get<std::string>("subgrid set name", "surface");

    std::size_t pos = mesh_plist.name().find('_');
    std::string parent_domain_name = mesh_plist.name().substr(pos+1,mesh_plist.name().size());
    auto parent = S.GetMesh(parent_domain_name);
    auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshSurfaceCell(parent, surface_setname));

    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(mesh_name, mesh, deformable);

  } else if (mesh_type == "domain set") {
    if (Amanzi::Keys::ends_with(mesh_name, "_*")) mesh_name = mesh_name.substr(0,mesh_name.length()-2);
    
    Teuchos::ParameterList& ds_list = mesh_plist.sublist("domain set parameters");
    auto parent_mesh_name = ds_list.get<std::string>("parent domain");
    auto parent_mesh = S.GetMesh(parent_mesh_name);
    auto entity_sets = ds_list.get<Teuchos::Array<std::string>>("entity set names").toVector();
    auto entity_kind = Amanzi::AmanziMesh::entity_kind(ds_list.get<std::string>("entity kind"));
    auto ds = Teuchos::rcp(new Amanzi::AmanziMesh::DomainSet(parent_mesh, entity_sets,
            entity_kind, mesh_name));
    S.RegisterDomainSet(mesh_name, ds);

    // a flyweight allows each subgrid model to share (topologically and
    // geometrically) identical meshes.
    bool flyweight = ds_list.get<bool>("flyweight mesh", false);

    // for each id in the regions of the parent mesh on entity, create a
    // subgrid mesh on MPI_COMM_SELF
    auto comm_self = Amanzi::getCommSelf();
    
    for (auto name_id : *ds) {
      Teuchos::ParameterList subgrid_i_list;
      if (ds_list.isSublist(name_id.first)) {
        subgrid_i_list = ds_list.sublist(name_id.first);
      } else {
        subgrid_i_list = ds_list.sublist(mesh_name+"_*");
      }
      subgrid_i_list.setName(name_id.first);
      auto& subgrid_i_param_list =
          subgrid_i_list.sublist(subgrid_i_list.get<std::string>("mesh type")+" parameters");
      if (!subgrid_i_param_list.isParameter("entity LID"))
        subgrid_i_param_list.set("entity LID", name_id.second);
      if (!subgrid_i_param_list.isParameter("entity kind"))
        subgrid_i_param_list.set("entity kind", Amanzi::AmanziMesh::entity_kind_string(ds->kind));
      if (!subgrid_i_param_list.isParameter("parent domain"))
        subgrid_i_param_list.set("parent domain", parent_mesh_name);
      createMesh(subgrid_i_list, comm_self, gm, S);
    }

  } else if (mesh_type == "Sperry 1D column") {
    // auto mesh = Testing::plantMesh(comm, gm, true);
    // bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    // checkVerifyMesh(mesh_plist, mesh);
    // S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name()), mesh, deformable);
    
  } else {
    Errors::Message msg;
    msg << "ATS Mesh Factory: unknown \"mesh type\" parameter \"" << mesh_type
        << "\" in mesh \"" << Amanzi::Keys::cleanPListName(mesh_plist.name()) << "\".";
    Exceptions::amanzi_throw(msg);
  }
}


bool checkVerifyMesh(Teuchos::ParameterList& mesh_plist,
                     Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh)
{
  // mesh verification
  AMANZI_ASSERT(!mesh.is_null());
  bool verify = mesh_plist.get<bool>("verify mesh", false);
  if (verify) {
    int num_procs = mesh->get_comm()->NumProc();
    int rank = mesh->get_comm()->MyPID();
    
    if (rank == 0)
      std::cout << "Verifying mesh with Mesh Audit..." << std::endl;
    if (num_procs == 1) {
      Amanzi::MeshAudit mesh_auditor(mesh);
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
      Amanzi::MeshAudit mesh_auditor(mesh, ofs);
      int status = mesh_auditor.Verify();        // check the mesh
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
  }  // if verify
  return true;
}


void
createMeshes(Teuchos::ParameterList& global_list,
             const Amanzi::Comm_ptr_type& comm,
             const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
             Amanzi::State& S)
{
  Teuchos::RCP<Teuchos::Time> volmeshtime =
      Teuchos::TimeMonitor::getNewCounter("volume mesh creation");
  Teuchos::TimeMonitor timer(*volmeshtime);
  
  Teuchos::ParameterList& meshes_list = global_list.sublist("mesh");

  // always try to do the domain mesh first
  if (meshes_list.isSublist("domain")) {
    createMesh(meshes_list.sublist("domain"), comm, gm, S);
  }

  // always try to do the surface mesh second
  if (meshes_list.isSublist("surface")) {
    createMesh(meshes_list.sublist("surface"), comm, gm, S);
  }

  // now do the rest
  for (auto sublist : meshes_list) {
    if (sublist.first != "domain" && sublist.first != "surface" &&
        meshes_list.isSublist(sublist.first)) {
      createMesh(meshes_list.sublist(sublist.first), comm, gm, S);
    }
  }

  // FIXME --etc
  // this should be dealt with somewhere else, and more generally
  // generalize vis for columns
  if (global_list.isSublist("visualization columns")) {
    auto surface_mesh = S.GetMesh("surface");
    Teuchos::ParameterList& vis_ss_plist = global_list.sublist("visualization columns"); 
    int nc = surface_mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
   
    for (int c=0; c!=nc; ++c){
      int id = surface_mesh->cell_map(false).GID(c);
      std::stringstream name_ss;
      name_ss << "column_" << id;
      vis_ss_plist.set("file name base", "visdump_"+name_ss.str());         
      global_list.set("visualization " +name_ss.str(), vis_ss_plist);
    }  
    global_list.remove("visualization columns");
  }
  
  // generalize vis for surface columns
  if (global_list.isSublist("visualization surface cells")) {
    auto surface_mesh = S.GetMesh("surface");
    Teuchos::ParameterList& vis_sf_plist = global_list.sublist("visualization surface cells");
    int nc = surface_mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    for (int c=0; c!=nc; ++c){
      int id = surface_mesh->cell_map(false).GID(c);
      std::stringstream name_ss, name_sf;
      name_sf << "surface_column_" << id;
      vis_sf_plist.set("file name base", "visdump_"+name_sf.str());
      global_list.set("visualization " +name_sf.str(), vis_sf_plist);
    }
    global_list.remove("visualization surface cells");
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //generalize checkpoint files for columns
  if(global_list.isSublist("checkpoints") && global_list.sublist("mesh").isSublist("column")){
  Teuchos::ParameterList& checkpoint_plist = global_list.sublist("checkpoints");
    std::stringstream name_check;
    name_check << rank;
    if (global_list.isSublist("checkpoints"))
      checkpoint_plist.set("file name base", "checkpoint_"+name_check.str() + "_");
    else
      checkpoint_plist.set("file name base", "checkpoint");
    global_list.set("checkpoint " +name_check.str(), checkpoint_plist);
    global_list.remove("checkpoints");
        
  }
  
  
  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();
  

}
} // namespace ATS
