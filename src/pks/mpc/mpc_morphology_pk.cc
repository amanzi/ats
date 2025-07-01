/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  This is the mpc_pk component of the Amanzi code.

*/

#include "mpc_morphology_pk.hh"
#include "PK_Helpers.hh"
#include "Mesh.hh"

namespace Amanzi {

Morphology_PK::Morphology_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                             const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree_or_fe_list, global_list, S, soln),
    MPCFlowTransport(pk_tree_or_fe_list, global_list, S, soln)
{
  // Create verbosity object.
  vo_ = Teuchos::null;
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = plist_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("Morphology_PK", vlist));
  domain_ = plist_->get<std::string>("domain name", "domain");
  name_ = "morphology pk";

  Teuchos::Array<std::string> pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
}


void
Morphology_PK::set_tags(const Tag& current, const Tag& next)
{

  MPCFlowTransport::set_tags(current, next);

}

void
Morphology_PK::parseParameterList(){

  MPCFlowTransport::parseParameterList();

  elevation_increase_key_ = Keys::getKey(domain_, "deformation");
  porosity_key_ = Keys::readKey(*plist_, domain_, "soil porosity", "soil_porosity");
  elev_key_ = Keys::readKey(*plist_, domain_, "elevation", "elevation");

  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  Teuchos::ParameterList& dens_list_next = S_->GetEvaluatorList(Keys::getKey("surface-mass_density_liquid", transport_next_tag));
  if (!dens_list_next.isParameter("evaluator type")) {
    dens_list_next.set<std::string>("evaluator type", "temporal interpolation");
    dens_list_next.set<std::string>("current tag", flow_current_tag.get());
    dens_list_next.set<std::string>("next tag", flow_next_tag.get());
  }

  Teuchos::ParameterList& sat_list_next = S_->GetEvaluatorList(Keys::getKey("surface-ponded_depth", transport_next_tag));
  if (!sat_list_next.isParameter("evaluator type")) {
    sat_list_next.set<std::string>("evaluator type", "temporal interpolation");
    sat_list_next.set<std::string>("current tag", flow_current_tag.get());
    sat_list_next.set<std::string>("next tag", flow_next_tag.get());
  }

  Teuchos::ParameterList& velo_list = S_->GetEvaluatorList(Keys::getKey("surface-velocity", transport_next_tag));
  if (!velo_list.isParameter("evaluator type")) {
    velo_list.set<std::string>("evaluator type", "alias");
    velo_list.set<std::string>("target", Keys::getKey("surface-velocity", flow_next_tag, true));
  }

  // Teuchos::ParameterList& elev_inc_list = S_->GetEvaluatorList(Keys::getKey(elevation_increase_key_, transport_next_tag));
  // if (!elev_inc_list.isParameter("evaluator type")) {
  //   elev_inc_list.set<std::string>("evaluator type", "alias");
  //   elev_inc_list.set<std::string>("target", Keys::getKey(elevation_increase_key_, flow_next_tag, true));
  // }

  dt_MPC_ = plist_->get<double>("dt MPC", 31557600);
  MSF_ = plist_->get<double>("morphological scaling factor", 1);

  mesh_ = S_->GetDeformableMesh(domain_);
  if (domain_ == "surface") {
    domain_3d_ = "surface_3d";
    domain_ss_ = "domain";

    vertex_coord_key_3d_ = Keys::getKey(domain_3d_, "vertex_coordinates");
    vertex_coord_key_ss_ = Keys::getKey(domain_ss_, "vertex_coordinates");

    surf3d_mesh_ = S_->GetDeformableMesh(domain_ + "_3d");
    mesh_ss_ = S_->GetDeformableMesh(domain_ss_);
  }

  if (!S_->FEList().isSublist(elev_key_)){
    S_->GetEvaluatorList(elev_key_).set("evaluator type", "meshed elevation");
  }
  S_->GetEvaluatorList(elev_key_).set("dynamic mesh", true);
  S_->GetEvaluatorList(elev_key_).set("deformation indicator", elevation_increase_key_);

}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
Morphology_PK::get_dt()
{
  if (dt_MPC_ < 0) {
    double dt = Amanzi::MPCFlowTransport::get_dt();
    set_dt(dt);
    return dt;
  } else {
    return dt_MPC_;
  }
}


void
Morphology_PK::Setup()
{

  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  Amanzi::MPCFlowTransport::Setup();

  // create storage for the vertex coordinates
  // we need to checkpoint those to be able to create
  // the deformed mesh after restart
  std::vector<AmanziMesh::Entity_kind> location(1);
  std::vector<int> num_dofs(1);
  std::vector<std::string> name(1);

  int dim = mesh_->getSpaceDimension();
  location[0] = AmanziMesh::NODE;
  num_dofs[0] = dim;
  name[0] = "node";

  S_->Require<CompositeVector, CompositeVectorSpace>(elevation_increase_key_, Tags::NEXT, elevation_increase_key_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAssign( elevation_increase_key_, Tags::NEXT, *S_, elevation_increase_key_);


  S_->Require<CompositeVector, CompositeVectorSpace>(porosity_key_, Tags::NEXT, porosity_key_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorAssign(porosity_key_, Tags::NEXT, *S_, porosity_key_);

  int num_veg_species = plist_->get<int>("number of vegitation species", 1);
  Key biomass_key = Keys::getKey(domain_, "biomass");
  Key stem_density_key = Keys::getKey(domain_, "stem_density");
  Key stem_height_key = Keys::getKey(domain_, "stem_height");
  Key stem_diameter_key = Keys::getKey(domain_, "stem_diameter");
  Key plant_area_key = Keys::getKey(domain_, "plant_area");

  location[0] = AmanziMesh::CELL;
  num_dofs[0] = num_veg_species;
  name[0] = "cell";

  requireEvaluatorAtNext(biomass_key, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponents(name, location, num_dofs);

  requireEvaluatorAtNext(plant_area_key, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponents(name, location, num_dofs);

  requireEvaluatorAtNext(stem_diameter_key, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponents(name, location, num_dofs);

  requireEvaluatorAtNext(stem_height_key, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponents(name, location, num_dofs);

  requireEvaluatorAtNext(stem_density_key, Tags::NEXT, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponents(name, location, num_dofs);

  S_->Require<double>("msl", Tags::NEXT);

}

 void
 Morphology_PK::Initialize()
{
  Amanzi::MPCFlowTransport::Initialize();

  // initialize the vertex coordinate of existing meshes

  // if (S_->HasRecord(vertex_coord_key_, Tags::NEXT)) Initialize_MeshVertices_(S_.ptr(), mesh_, vertex_coord_key_, Tags::NEXT);
  // if (S_->HasRecord(vertex_coord_key_3d_, Tags::NEXT)) Initialize_MeshVertices_(S_.ptr(), surf3d_mesh_, vertex_coord_key_3d_, Tags::NEXT);
  // if (S_->HasRecord(vertex_coord_key_ss_, Tags::NEXT)) Initialize_MeshVertices_(S_.ptr(), mesh_ss_, vertex_coord_key_ss_, Tags::NEXT);

  if (S_->HasRecord(elevation_increase_key_, Tags::NEXT)) {
    S_->GetW<CompositeVector>(elevation_increase_key_, Tags::NEXT, elevation_increase_key_).PutScalar(0.);
    S_->GetRecordW(elevation_increase_key_, Tags::NEXT, elevation_increase_key_).set_initialized();
  }

  flow_pk_ = Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[0]);
  sed_transport_pk_ = sub_pks_[1];


  const Epetra_MultiVector& dz =
    *S_->Get<CompositeVector>(elevation_increase_key_, Tags::NEXT).ViewComponent("cell", false);

  dz_accumul_ = Teuchos::rcp(new Epetra_MultiVector(dz));
  dz_accumul_->PutScalar(0.);
}

void
Morphology_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  Amanzi::MPCFlowTransport::CommitStep(t_old, t_new, tag);

  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  Key slope_key = Keys::readKey(*plist_, domain_, "slope magnitude", "slope_magnitude");

  bool chg = S_->GetEvaluator(elev_key, Tags::NEXT).Update(*S_, name_);


  Key biomass_key = Keys::getKey(domain_, "biomass");
  chg = S_->GetEvaluator(biomass_key, Tags::NEXT).Update(*S_, name_);
}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
Morphology_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  Teuchos::OSTab tab = vo_->getOSTab();

  //Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  Epetra_MultiVector& dz =
    *S_->GetW<CompositeVector>(elevation_increase_key_, Tags::NEXT, elevation_increase_key_).ViewComponent("cell", false);
  dz.PutScalar(0.);

  fail = Amanzi::MPCFlowTransport::AdvanceStep(t_old, t_new, reinit);

  if (!fail) {

    double max_dz, min_dz;
    dz.MinValue(&min_dz);
    dz.MaxValue(&max_dz);

    dz.Scale(MSF_);

    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Deformation min " << min_dz << " max " << max_dz << "\n";

    //    dz_accumul_->Update(1, dz, 1);
    Update_MeshVertices_(S_.ptr(), Tags::NEXT);

  }

  return fail;

}

void
Morphology_PK::Initialize_MeshVertices_(const Teuchos::Ptr<State>& S,
                                        Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                        Key vert_field_key, Tag field_tag)
{
  // spatial dimension
  int dim = mesh->getSpaceDimension();
  Amanzi::AmanziGeometry::Point coords(dim);
  // number of vertices
  int nV = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE, Amanzi::AmanziMesh::Parallel_kind::OWNED);

  Epetra_MultiVector& vc =
    *S->GetPtrW<CompositeVector>(vert_field_key, field_tag, "state")->ViewComponent("node", false);

  // search the id of the mid point on the top
  for (int iV = 0; iV < nV; iV++) {
    // get the coords of the node
    coords = mesh->getNodeCoordinate(iV);
    for (int s = 0; s < dim; ++s) { vc[s][iV] = coords[s]; }
  }

  S->GetPtrW<CompositeVector>(vert_field_key, field_tag, "state")->ScatterMasterToGhosted("node");
  S->GetRecordW(vert_field_key, field_tag, "state").set_initialized();
}

void
Morphology_PK::Update_MeshVertices_(const Teuchos::Ptr<State>& S, const Tag& tag)
{
  // spatial dimension
  int dim = mesh_ss_->getSpaceDimension();
  Amanzi::AmanziGeometry::Point coords(dim);
  // number of vertices

  Epetra_MultiVector& vc =
    *S->GetPtrW<CompositeVector>(vertex_coord_key_ss_, tag, vertex_coord_key_ss_)->ViewComponent("node", true);

  Epetra_MultiVector& vc_surf_3d =
    *S->GetPtrW<CompositeVector>(vertex_coord_key_3d_, tag, vertex_coord_key_3d_)->ViewComponent("node", true);


  const Epetra_MultiVector& dz =
    *S->Get<CompositeVector>(elevation_increase_key_, tag).ViewComponent("cell");

  int nsurf_cells = dz.MyLength();
  int nsurfnodes = surf3d_mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE,
                                                       AmanziMesh::Parallel_kind::ALL);
  int nnodes = mesh_ss_->getNumEntities(AmanziMesh::Entity_kind::NODE,
                                                       AmanziMesh::Parallel_kind::ALL);

  AmanziMesh::Entity_ID_View surface_nodeids("surface_nodeids", nsurfnodes);
  AmanziMesh::Point_View surface_newpos("surface_newpos", nsurfnodes);
  AmanziMesh::Entity_ID_View nodeids("nodeids", nnodes);
  AmanziMesh::Point_View newpos("newpos", nnodes);

  for (int c = 0; c < nsurf_cells; c++) {
    AmanziMesh::Entity_ID domain_face;
    domain_face = mesh_->getEntityParent(AmanziMesh::CELL, c);

    auto nodes = mesh_ss_->getFaceNodes(domain_face);
    int nnodes = nodes.size();
    for (int i = 0; i < nnodes; i++) {
      mesh_ss_->getNodeCoordinate(nodes[i]);
      mesh_ss_->getNodeCells(nodes[i]);
      auto cells = mesh_ss_->getNodeCells(nodes[i]);
      int nsize = cells.size();

      vc[2][nodes[i]] += dz[0][c] / nsize;
    }
  }

  S_->GetPtrW<CompositeVector>(vertex_coord_key_ss_, tag, vertex_coord_key_ss_)->ScatterMasterToGhosted("node");

  for (int i=0; i<nnodes; ++i){
    nodeids[i] = i;
    coords = mesh_ss_->getNodeCoordinate(i);
    coords[2] = vc[2][i];
    newpos[i] = coords;
  }

  AmanziMesh::deform(*mesh_ss_, nodeids, newpos);


  for (int c = 0; c < nsurf_cells; c++) {
    AmanziMesh::Entity_ID domain_face;

    auto nodes = surf3d_mesh_->getFaceNodes(c);
    int nnodes = nodes.size();
    for (int i = 0; i < nnodes; i++) {
      coords = surf3d_mesh_->getNodeCoordinate(nodes[i]);
      auto cells = surf3d_mesh_->getNodeCells(nodes[i]);
      int nsize = cells.size();
      //coords[2] += dz / nsize;
      vc_surf_3d[2][nodes[i]] += dz[0][c]/ nsize;

    }
  }

  S_->GetPtrW<CompositeVector>(vertex_coord_key_3d_, tag, vertex_coord_key_3d_)->ScatterMasterToGhosted("node");

  for (int i=0; i<nsurfnodes; ++i){
    surface_nodeids[i] = i;
    coords = surf3d_mesh_->getNodeCoordinate(i);
    coords[2] = vc_surf_3d[2][i];
    surface_newpos[i] = coords;
  }

  // AmanziMesh::deform(*surf3d_mesh_, surface_nodeids, surface_newpos);

  deform_eval_ =
     Teuchos::rcp_dynamic_cast<EvaluatorPrimaryCV>(S_->GetEvaluatorPtr(elevation_increase_key_, tag));
  deform_eval_->SetChanged();


}


} // namespace Amanzi
