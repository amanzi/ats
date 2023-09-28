/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A coupler which integrates surface and subsurface flow implicitly.
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "pk_helpers.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "mpc_coupled_water.hh"

namespace Amanzi {

MPCCoupledWater::MPCCoupledWater(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_plist, S, soln),
    StrongMPC<PK_PhysicalBDF_Default>(pk_tree, global_plist, S, soln)
{}


void
MPCCoupledWater::Setup()
{
  // tweak the sub-PK parameter lists
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");

  // -- turn on coupling
  pks_list_->sublist(names[0]).set("coupled to surface via flux", true);
  pks_list_->sublist(names[1]).set("coupled to subsurface via flux", true);

  // -- ensure local ops are suface ops
  pks_list_->sublist(names[1]).sublist("diffusion preconditioner").set("surface operator", true);
  pks_list_->sublist(names[1]).sublist("accumulation preconditioner").set("surface operator", true);

  domain_ss_ = pks_list_->sublist(names[0]).get<std::string>("domain name", "domain");
  domain_surf_ = pks_list_->sublist(names[1]).get<std::string>("domain name", "surface");

  // grab the meshes
  surf_mesh_ = S_->GetMesh(domain_surf_);
  domain_mesh_ = S_->GetMesh(domain_ss_);

  // cast the PKs
  domain_flow_pk_ = sub_pks_[0];
  surf_flow_pk_ = sub_pks_[1];

  // call the MPC's setup, which calls the sub-pk's setups
  StrongMPC<PK_PhysicalBDF_Default>::Setup();

  // require the coupling fields, claim ownership
  exfilt_key_ =
    Keys::readKey(*plist_, domain_surf_, "exfiltration flux", "surface_subsurface_flux");
  S_->Require<CompositeVector, CompositeVectorSpace>(exfilt_key_, tag_next_, exfilt_key_)
    .SetMesh(surf_mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  requireEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // Create the preconditioner.
  // -- collect the preconditioners
  precon_ = domain_flow_pk_->preconditioner();
  precon_surf_ = surf_flow_pk_->preconditioner();

  // -- set parameters for an inverse
  Teuchos::ParameterList inv_list = plist_->sublist("inverse");
  inv_list.setParameters(plist_->sublist("preconditioner"));
  inv_list.setParameters(plist_->sublist("linear solver"));
  precon_->set_inverse_parameters(inv_list);

  // -- push the surface local ops into the subsurface global operator
  for (Operators::Operator::op_iterator op = precon_surf_->begin(); op != precon_surf_->end();
       ++op) {
    precon_->OpPushBack(*op);
  }

  // set up the Water delegate
  Teuchos::RCP<Teuchos::ParameterList> water_list = Teuchos::sublist(plist_, "water delegate");
  water_ = Teuchos::rcp(new MPCDelegateWater(water_list, S_, domain_ss_, domain_surf_));
  water_->set_tags(tag_current_, tag_next_);
  water_->set_indices(0, 1);

  // grab the debuggers
  domain_db_ = domain_flow_pk_->debugger();
  water_->set_db(domain_db_);
  surf_db_ = surf_flow_pk_->debugger();

  // With this MPC, thanks to the form of the error/solver, it is often easier
  // to figure out what subsurface face or cell to debug, and it is hard to
  // figure out what the corresponding cell of the surface system is.
  // Therefore, for all debug cells of the subsurface, if that cell is in the
  // top layer of cells, we add the corresponding face's surface cell.
  AmanziMesh::Entity_ID_List debug_cells = domain_db_->get_cells();

  int ncells_surf =
    surf_mesh_->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  if (debug_cells.size() > 0) {
    const auto& domain_cell_map = domain_mesh_->cell_map(false);
    const auto& surf_cell_map = surf_mesh_->cell_map(false);
    AmanziMesh::Entity_ID_List surf_debug_cells;
    for (int sc = 0; sc != ncells_surf; ++sc) {
      int f = surf_mesh_->entity_get_parent(AmanziMesh::Entity_kind::CELL, sc);
      AmanziMesh::Entity_ID_List fcells;
      domain_mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &fcells);
      AMANZI_ASSERT(fcells.size() == 1);
      auto gid = domain_cell_map.GID(fcells[0]);
      if (std::find(debug_cells.begin(), debug_cells.end(), gid) != debug_cells.end())
        surf_debug_cells.emplace_back(surf_cell_map.GID(sc));
    }
    if (surf_debug_cells.size() > 0) surf_db_->add_cells(surf_debug_cells);
  }
}


void
MPCCoupledWater::Initialize()
{
  // initialize coupling terms
  S_->GetPtrW<CompositeVector>(exfilt_key_, tag_next_, exfilt_key_)->PutScalar(0.);
  S_->GetRecordW(exfilt_key_, tag_next_, exfilt_key_).set_initialized();
  changedEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // Initialize all sub PKs.
  MPC<PK_PhysicalBDF_Default>::Initialize();

  // ensure continuity of ICs... subsurface takes precedence.
  CopySubsurfaceToSurface(S_->Get<CompositeVector>(Keys::getKey(domain_ss_, "pressure"), tag_next_),
                          S_->GetW<CompositeVector>(Keys::getKey(domain_surf_, "pressure"),
                                                    tag_next_,
                                                    sub_pks_[1]->name()));

  // Initialize my timestepper.
  PK_BDF_Default::Initialize();
}


// -- computes the non-linear functional g = g(t,u,udot)
//    By default this just calls each sub pk FunctionalResidual().
void
MPCCoupledWater::FunctionalResidual(double t_old,
                                    double t_new,
                                    Teuchos::RCP<TreeVector> u_old,
                                    Teuchos::RCP<TreeVector> u_new,
                                    Teuchos::RCP<TreeVector> g)
{
  // propagate updated info into state
  Solution_to_State(*u_new, tag_next_);

  // Evaluate the surface flow residual
  surf_flow_pk_->FunctionalResidual(
    t_old, t_new, u_old->SubVector(1), u_new->SubVector(1), g->SubVector(1));

  // The residual of the surface flow equation provides the water flux from
  // subsurface to surface.
  Epetra_MultiVector& source =
    *S_->GetW<CompositeVector>(exfilt_key_, tag_next_, exfilt_key_).ViewComponent("cell", false);
  source = *g->SubVector(1)->Data()->ViewComponent("cell", false);
  changedEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->FunctionalResidual(
    t_old, t_new, u_old->SubVector(0), u_new->SubVector(0), g->SubVector(0));

  // All surface to subsurface fluxes have been taken by the subsurface.
  g->SubVector(1)->Data()->ViewComponent("cell", false)->PutScalar(0.);
}


// -- Apply preconditioner to u and returns the result in Pu.
int
MPCCoupledWater::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Precon application:" << std::endl;

  // call the precon's inverse
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying subsurface operator." << std::endl;
  int ierr = precon_->ApplyInverse(*u->SubVector(0)->Data(), *Pu->SubVector(0)->Data());

  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying  CopySubsurfaceToSurface." << std::endl;
  // Copy subsurface face corrections to surface cell corrections
  CopySubsurfaceToSurface(*Pu->SubVector(0)->Data(), *Pu->SubVector(1)->Data());

  // // Derive surface face corrections.
  // UpdateConsistentFaceCorrectionWater_(u, Pu);

  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");

    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(u->SubVector(0)->Data().ptr());
    vecs.push_back(Pu->SubVector(0)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = u->SubVector(1)->Data().ptr();
    vecs[1] = Pu->SubVector(1)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);
  }
  return (ierr > 0) ? 0 : 1;
}


// -- Modify the predictor.
bool
MPCCoupledWater::ModifyPredictor(double h,
                                 Teuchos::RCP<const TreeVector> u0,
                                 Teuchos::RCP<TreeVector> u)
{
  bool modified = false;

  // Calculate consistent faces
  modified = domain_flow_pk_->ModifyPredictor(h, u0->SubVector(0), u->SubVector(0));

  // Merge surface cells with subsurface faces
  if (modified) {
    Teuchos::RCP<const CompositeVector> h_prev =
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_surf_, "ponded_depth"), tag_next_);
    MergeSubsurfaceAndSurfacePressure(*h_prev, *u->SubVector(0)->Data(), *u->SubVector(1)->Data());
  }

  // Hack surface faces
  bool newly_modified = false;
  bool has_face = u->SubVector(0)->Data()->HasComponent("face");
  if (has_face) {
    newly_modified |= water_->ModifyPredictor_Heuristic<AmanziMesh::FACE>(h, u);
    newly_modified |= water_->ModifyPredictor_WaterSpurtDamp<AmanziMesh::FACE>(h, u);
  } else {
    newly_modified |= water_->ModifyPredictor_Heuristic<AmanziMesh::BOUNDARY_FACE>(h, u);
    newly_modified |= water_->ModifyPredictor_WaterSpurtDamp<AmanziMesh::BOUNDARY_FACE>(h, u);
  }
  modified |= newly_modified;

  // -- copy surf --> sub
  if (newly_modified) {
    CopySurfaceToSubsurface(*u->SubVector(1)->Data(), *u->SubVector(0)->Data());
  }

  // Calculate consistent surface faces
  modified |= surf_flow_pk_->ModifyPredictor(h, u0->SubVector(1), u->SubVector(1));
  return modified;
}


// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCCoupledWater::ModifyCorrection(double h,
                                  Teuchos::RCP<const TreeVector> res,
                                  Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "NKA'd, PC'd correction." << std::endl;

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");

    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(res->SubVector(0)->Data().ptr());
    vecs.push_back(du->SubVector(0)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = res->SubVector(1)->Data().ptr();
    vecs[1] = du->SubVector(1)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);
  }

  // modify correction using sub-pk approaches
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult modified_res =
    StrongMPC<PK_PhysicalBDF_Default>::ModifyCorrection(h, res, u, du);

  // modify correction using water approaches
  int n_modified = 0;
  double damping;
  bool has_face = u->SubVector(0)->Data()->HasComponent("face");
  if (has_face) {
    n_modified += water_->ModifyCorrection_WaterFaceLimiter<AmanziMesh::FACE>(h, res, u, du);
    double damping1 = water_->ModifyCorrection_WaterSpurtDamp<AmanziMesh::FACE>(h, res, u, du);
    double damping2 =
      water_->ModifyCorrection_DesaturatedSpurtDamp<AmanziMesh::FACE>(h, res, u, du);
    damping = std::min(damping1, damping2);
    n_modified += water_->ModifyCorrection_WaterSpurtCap<AmanziMesh::FACE>(h, res, u, du, damping);
    n_modified +=
      water_->ModifyCorrection_DesaturatedSpurtCap<AmanziMesh::FACE>(h, res, u, du, damping);
  } else {
    n_modified +=
      water_->ModifyCorrection_WaterFaceLimiter<AmanziMesh::BOUNDARY_FACE>(h, res, u, du);
    double damping1 =
      water_->ModifyCorrection_WaterSpurtDamp<AmanziMesh::BOUNDARY_FACE>(h, res, u, du);
    double damping2 =
      water_->ModifyCorrection_DesaturatedSpurtDamp<AmanziMesh::BOUNDARY_FACE>(h, res, u, du);
    damping = std::min(damping1, damping2);
    n_modified +=
      water_->ModifyCorrection_WaterSpurtCap<AmanziMesh::BOUNDARY_FACE>(h, res, u, du, damping);
    n_modified += water_->ModifyCorrection_DesaturatedSpurtCap<AmanziMesh::BOUNDARY_FACE>(
      h, res, u, du, damping);
  }

  // -- accumulate globally
  int n_modified_l = n_modified;
  u->SubVector(0)->Data()->Comm()->SumAll(&n_modified_l, &n_modified, 1);
  bool modified = (n_modified > 0) || (damping < 1.);

  // -- calculate consistent subsurface cells
  if (modified) {
    // if (consistent_cells_) {
    //   // Derive subsurface cell corrections.
    //   precon_->UpdateConsistentCellCorrection(
    //       *u->SubVector(0)->Data(),
    //       du->SubVector(0)->Data().ptr());
    // }

    // Copy subsurface face corrections to surface cell corrections
    CopySubsurfaceToSurface(*du->SubVector(0)->Data(), *du->SubVector(1)->Data());
  }

  // if (modified) {
  //   // Derive surface face corrections.
  //   UpdateConsistentFaceCorrectionWater_(res, du);
  // }

  // dump to screen
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Modified correction." << std::endl;

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");

    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(res->SubVector(0)->Data().ptr());
    vecs.push_back(du->SubVector(0)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = res->SubVector(1)->Data().ptr();
    vecs[1] = du->SubVector(1)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);
  }

  return (modified_res || modified) ?
           AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING :
           AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}


double
MPCCoupledWater::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> res)
{
  // move the surface face residual onto the surface cell.
  auto res2 = Teuchos::rcp(new TreeVector(*res, INIT_MODE_COPY));
  auto& res_face = *res2->SubVector(0)->Data()->ViewComponent("face", false);
  auto& res_surf_cell = *res2->SubVector(1)->Data()->ViewComponent("cell", false);
  const auto& u_surf_cell = *u->SubVector(1)->Data()->ViewComponent("cell", false);
  double p_atm = S_->Get<double>("atmospheric_pressure", Tags::NEXT);
  for (int c = 0; c != u_surf_cell.MyLength(); ++c) {
    if (u_surf_cell[0][c] > p_atm) {
      auto f = surf_mesh_->entity_get_parent(AmanziMesh::Entity_kind::CELL, c);
      res_surf_cell[0][c] = res_face[0][f];
      res_face[0][f] = 0.;
    }
  }
  return StrongMPC<PK_PhysicalBDF_Default>::ErrorNorm(u, res2);
}


// void
// MPCCoupledWater::UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
//         const Teuchos::RCP<TreeVector>& Pu) {
//   Teuchos::OSTab tab = vo_->getOSTab();

//   Teuchos::RCP<CompositeVector> surf_Pp = Pu->SubVector(1)->Data();
//   Epetra_MultiVector& surf_Pp_c = *surf_Pp->ViewComponent("cell",false);

//   Teuchos::RCP<const CompositeVector> surf_p = u->SubVector(1)->Data();

//   // Calculate delta h on the surface
//   Teuchos::RCP<CompositeVector> surf_Ph = Teuchos::rcp(new CompositeVector(*surf_Pp));
//   surf_Ph->PutScalar(0.);

//   // old ponded depth
//   S_next_->GetEvaluator("ponded_depth")->Update(S_next_.ptr(), name_);
//   *surf_Ph->ViewComponent("cell",false) = *S_next_->Get<CompositeVector>("ponded_depth").ViewComponent("cell",false);

//   // new ponded depth
//   Teuchos::RCP<TreeVector> tv_p = Teuchos::rcp(new TreeVector());
//   Teuchos::RCP<CompositeVector> cv_p = S_next_->GetPtrW<CompositeVector>("surface-pressure", sub_pks_[1]->name());
//   cv_p->ViewComponent("cell",false)->Update(-1., surf_Pp_c, 1.);
//   tv_p->SetData(cv_p);

//   sub_pks_[1]->ChangedSolution();

//   if (sub_pks_[1]->IsAdmissible(tv_p)) {
//     S_next_->GetEvaluator("ponded_depth")->Update(S_next_.ptr(), name_);

//     // put delta ponded depth into surf_Ph_cell
//     surf_Ph->ViewComponent("cell",false)
//         ->Update(-1., *S_next_->Get<CompositeVector>("ponded_depth").ViewComponent("cell",false), 1.);

//     // update delta faces
//     precon_surf_->UpdateConsistentFaceCorrection(*surf_p, surf_Ph.ptr());
//     *surf_Pp->ViewComponent("face",false) = *surf_Ph->ViewComponent("face",false);

//     // dump to screen
//     if (vo_->os_OK(Teuchos::VERB_HIGH)) {
//       *vo_->os() << "Precon water correction." << std::endl;

//       std::vector<std::string> vnames;
//       vnames.push_back("pd_new");
//       vnames.push_back("delta_pd");

//       std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
//       vecs.push_back(S_next_->GetPtr<CompositeVector>("ponded_depth").ptr());
//       vecs.push_back(surf_Ph.ptr());
//       surf_db_->WriteVectors(vnames, vecs, true);
//     }
//   }

//   // revert solution so we don't break things
//   S_next_->GetPtrW<CompositeVector>("surface-pressure",sub_pks_[1]->name())
//       ->ViewComponent("cell",false)->Update(1., surf_Pp_c, 1.);
//   sub_pks_[1]->ChangedSolution();
// }


} // namespace Amanzi
