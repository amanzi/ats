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
{}


void
Morphology_PK::parseParameterList()
{
  MPCFlowTransport::parseParameterList();

  elevation_increase_key_ = Keys::readKey(*plist_, domain_, "deformation", "deformation");
  porosity_key_ = Keys::readKey(*plist_, domain_, "soil porosity", "soil_porosity");
  elev_key_ = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  dens_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
  velo_key_ = Keys::readKey(*plist_, domain_, "velocity", "velocity");

  biomass_key_ = Keys::readKey(*plist_, domain_, "biomass", "biomass");
  stem_density_key_ = Keys::readKey(*plist_, domain_, "stem density", "stem_density");
  stem_height_key_ = Keys::readKey(*plist_, domain_, "stem height", "stem_height");
  stem_diameter_key_ = Keys::readKey(*plist_, domain_, "stem diameter", "stem_diameter");
  plant_area_key_ = Keys::readKey(*plist_, domain_, "plant area", "plant_area");

  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  Teuchos::ParameterList& dens_list_next = S_->GetEvaluatorList(Keys::getKey(dens_key_, transport_next_tag));
  if (!dens_list_next.isParameter("evaluator type")) {
    dens_list_next.set<std::string>("evaluator type", "temporal interpolation");
    dens_list_next.set<std::string>("current tag", flow_current_tag.get());
    dens_list_next.set<std::string>("next tag", flow_next_tag.get());
  }

  Teuchos::ParameterList& velo_list = S_->GetEvaluatorList(Keys::getKey(velo_key_, transport_next_tag));
  if (!velo_list.isParameter("evaluator type")) {
    velo_list.set<std::string>("evaluator type", "alias");
    velo_list.set<std::string>("target", Keys::getKey(velo_key_, flow_next_tag, true));
  }

  MSF_ = plist_->get<double>("morphological scaling factor", 1);

  mesh_ = S_->GetDeformableMesh(domain_);

  domain_3d_ = Keys::readDomainHint(*plist_, domain_, "surface", "surface_3d");
  domain_ss_ = Keys::readDomainHint(*plist_, domain_, "surface", "domain");

  vertex_coord_key_3d_ = Keys::readKey(*plist_, domain_3d_, "vertex coordinates", "vertex_coordinates");
  vertex_coord_key_ss_ = Keys::readKey(*plist_, domain_ss_, "vertex coordinates", "vertex_coordinates");

  surf3d_mesh_ = S_->GetDeformableMesh(domain_3d_);
  mesh_ss_ = S_->GetDeformableMesh(domain_ss_);

  if (!S_->HasEvaluatorList(elev_key_)) {
    Teuchos::ParameterList& elev_list = S_->GetEvaluatorList(elev_key_);
    elev_list.set("evaluator type", "meshed elevation");
    elev_list.set("dynamic mesh", true);
    elev_list.set("deformation indicator", elevation_increase_key_);
  }
}


void
Morphology_PK::Setup()
{
  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  Amanzi::MPCFlowTransport::Setup();

  requireEvaluatorAssign( elevation_increase_key_, tag_next_, *S_, elevation_increase_key_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireEvaluatorAssign(porosity_key_, tag_next_, *S_, porosity_key_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  int num_veg_species = plist_->get<int>("number of vegetation species", 1);

  requireEvaluatorAtNext(biomass_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_veg_species);

  requireEvaluatorAtNext(plant_area_key, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_veg_species);

  requireEvaluatorAtNext(stem_diameter_key, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_veg_species);

  requireEvaluatorAtNext(stem_height_key, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_veg_species);

  requireEvaluatorAtNext(stem_density_key, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_veg_species);

  S_->Require<double>("msl", tag_next_);
}

 void
 Morphology_PK::Initialize()
{
  Amanzi::MPCFlowTransport::Initialize();

  if (S_->HasRecord(elevation_increase_key_, tag_next_)) {
    S_->GetW<CompositeVector>(elevation_increase_key_, tag_next_, elevation_increase_key_).PutScalar(0.);
    S_->GetRecordW(elevation_increase_key_, tag_next_, elevation_increase_key_).set_initialized();
  }

  flow_pk_ = Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[0]);
  sed_transport_pk_ = sub_pks_[1];


  const Epetra_MultiVector& dz =
    *S_->Get<CompositeVector>(elevation_increase_key_, tag_next_).ViewComponent("cell", false);

  dz_accumul_ = Teuchos::rcp(new Epetra_MultiVector(dz));
  dz_accumul_->PutScalar(0.);
}


void
Morphology_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  Amanzi::MPCFlowTransport::CommitStep(t_old, t_new, tag);

  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  Key slope_key = Keys::readKey(*plist_, domain_, "slope magnitude", "slope_magnitude");

  bool chg = S_->GetEvaluator(elev_key, tag_next_).Update(*S_, name_);


  Key biomass_key = Keys::getKey(domain_, "biomass");
  chg = S_->GetEvaluator(biomass_key, tag_next_).Update(*S_, name_);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
Morphology_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  Teuchos::OSTab tab = vo_->getOSTab();

  S_->GetW<CompositeVector>(elevation_increase_key_, tag_next_, elevation_increase_key_).PutScalar(0.);
  fail = Amanzi::MPCFlowTransport::AdvanceStep(t_old, t_new, reinit);

  if (!fail) {
    const Epetra_MultiVector& dz = S_->Get<CompositeVector>(elevation_increase_key_, tag_next_, elevation_increase_key_)
      .ViewComponent("cell", false);
    double max_dz, min_dz;
    dz.MinValue(&min_dz);
    dz.MaxValue(&max_dz);

    dz.Scale(MSF_);

    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Deformation min " << min_dz << " max " << max_dz << "\n";

    Update_MeshVertices_(S_.ptr(), tag_next_);
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
