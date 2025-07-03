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

  domain_ = plist_->get<std::string>("domain name");
  elevation_increase_key_ = Keys::readKey(*plist_, domain_, "deformation", "deformation");

  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  if (!S_->HasEvaluatorList(elev_key)) {
    Teuchos::ParameterList& elev_list = S_->GetEvaluatorList(elev_key);
    elev_list.set("evaluator type", "meshed elevation");
    elev_list.set("dynamic mesh", true);
    elev_list.set("deformation indicator", elevation_increase_key_);
  }

  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  dens_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
  pd_key_ = Keys::readKey(*plist_, domain_, "ponded depth", "ponded_depth");
  // if subcycling transport, need an interpolation eval for density
  if (flow_next_tag != transport_next_tag) {
    // flow's current is a copy
    requireEvaluatorAtCurrent(dens_key_, flow_current_tag, *S_, name_);
    requireEvaluatorAtCurrent(pd_key_, flow_current_tag, *S_, name_);

    // transport's next is an interpolation
    Teuchos::ParameterList& dens_list_next = S_->GetEvaluatorList(Keys::getKey(dens_key_, transport_next_tag));
    if (!dens_list_next.isParameter("evaluator type")) {
      dens_list_next.set<std::string>("evaluator type", "temporal interpolation");
      dens_list_next.set<std::string>("current tag", flow_current_tag.get());
      dens_list_next.set<std::string>("next tag", flow_next_tag.get());
    }
    Teuchos::ParameterList& pd_list_next = S_->GetEvaluatorList(Keys::getKey(pd_key_, transport_next_tag));
    if (!pd_list_next.isParameter("evaluator type")) {
      pd_list_next.set<std::string>("evaluator type", "temporal interpolation");
      pd_list_next.set<std::string>("current tag", flow_current_tag.get());
      pd_list_next.set<std::string>("next tag", flow_next_tag.get());
    }
  }

  MSF_ = plist_->get<double>("morphological scaling factor", 1);

  domain_3d_ = Keys::readDomainHint(*plist_, domain_, "surface", "surface_3d");
  domain_ss_ = Keys::readDomainHint(*plist_, domain_, "surface", "domain");

  vertex_coord_key_ = Keys::readKey(*plist_, domain_3d_, "vertex coordinates", "vertex_coordinates");

  mesh_ = S_->GetDeformableMesh(domain_);
  mesh_3d_ = S_->GetDeformableMesh(domain_3d_);
  mesh_ss_ = S_->GetDeformableMesh(domain_ss_);
}


void
Morphology_PK::Setup()
{
  auto [flow_current_tag, flow_next_tag] = tags_[0];
  auto [transport_current_tag, transport_next_tag] = tags_[1];

  Amanzi::MPCFlowTransport::Setup();

  requireEvaluatorAtNext(elevation_increase_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<CompositeVector, CompositeVectorSpace>(vertex_coord_key_, tag_next_, vertex_coord_key_)
    .SetMesh(mesh_3d_)
    ->SetGhosted()
    ->SetComponent("node", AmanziMesh::Entity_kind::NODE, mesh_3d_->getSpaceDimension());
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
Morphology_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  Teuchos::OSTab tab = vo_->getOSTab();

  fail = Amanzi::MPCFlowTransport::AdvanceStep(t_old, t_new, reinit);

  if (!fail) {
    const Epetra_MultiVector& dz = *S_->Get<CompositeVector>(elevation_increase_key_, tag_next_)
      .ViewComponent("cell", false);
    double max_dz, min_dz;
    dz.MinValue(&min_dz);
    dz.MaxValue(&max_dz);

    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Deformation min " << min_dz << " max " << max_dz << "\n";

    Update_MeshVertices_(S_.ptr(), tag_next_);
  }
  return fail;
}

void
Morphology_PK::Update_MeshVertices_(const Teuchos::Ptr<State>& S, const Tag& tag)
{
  S->Get<CompositeVector>(elevation_increase_key_, tag).ScatterMasterToGhosted("cell");
  { // This is done in a context to avoid problems with the scatter and views
    const Epetra_MultiVector& dz =
      *S->Get<CompositeVector>(elevation_increase_key_, tag).ViewComponent("cell", true);

    // NOTE: this will get much cleaner in Tpetra, where the view and the
    // vector can share the same memory!  For now we need the view for the
    // deform() call, and the vector for the scatter.
    //
    // NOTE that this only needs to exist on the surface mesh, but needs to be
    // the 3d coordinates.
    Epetra_MultiVector& vc =
      *S->GetW<CompositeVector>(vertex_coord_key_, tag, vertex_coord_key_).ViewComponent("node", false);

    for (int c = 0; c != dz.MyLength(); ++c) {
      auto surf_nodes = mesh_->getCellNodes(c);
      for (auto& n : surf_nodes) {
        if (n < vc.MyLength()) {
          int ncells = mesh_->getNodeCells(n, AmanziMesh::Parallel_kind::ALL).size();
          vc[2][n] += MSF_ * dz[0][c] / ncells;
        }
      }
    }
  }

  S_->Get<CompositeVector>(vertex_coord_key_, tag).ScatterMasterToGhosted("node");

  { // now that we have scattered, move the coordinates to a view and deform
    //
    // NOTE: must deform the ghost nodes as well as owned nodes!
    const Epetra_MultiVector& vc =
      *S->Get<CompositeVector>(vertex_coord_key_, tag).ViewComponent("node", true);

    // NOTE: nodeids will be subsurface mesh node ids!
    AmanziMesh::Entity_ID_View nodeids("nodeids", vc.MyLength());
    AmanziMesh::Entity_ID_View nodeids_surf("nodeids_surf", vc.MyLength());
    AmanziMesh::Mesh::Point_View newpos("newpos", vc.MyLength());

    AmanziMesh::cEntity_ID_View parent_ids = mesh_->getEntityParents(AmanziMesh::Entity_kind::NODE);
    for (int n = 0; n != vc.MyLength(); ++n) {
      AmanziMesh::Entity_ID ss_node = parent_ids[n];
      nodeids_surf[n] = n;
      nodeids[n] = ss_node;
      AmanziGeometry::Point p(vc[0][n], vc[1][n], vc[2][n]);
      newpos[n] = p;
    }

    AmanziMesh::deform(*mesh_ss_, nodeids, newpos);
    AmanziMesh::deform(*mesh_3d_, nodeids_surf, newpos);

    // zero out the dz field to reset for the next step
    S_->GetW<CompositeVector>(elevation_increase_key_, tag_next_,
            sub_pks_[1]->name()).PutScalar(0.);

    // mark the indicator field as changed now to force recomputation of
    // dependencies on deformed mesh
    changedEvaluatorPrimary(elevation_increase_key_, tag_next_, *S_);
  }
}


void
Morphology_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  MPCFlowTransport::CommitStep(t_old, t_new, tag);

  // also save our flow quantities, needed for interpolation
  auto [flow_current_tag, flow_next_tag] = tags_[0];
  assign(dens_key_, flow_current_tag, flow_next_tag, *S_);
  assign(pd_key_, flow_current_tag, flow_next_tag, *S_);
}



} // namespace Amanzi
