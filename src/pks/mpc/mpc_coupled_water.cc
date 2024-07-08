/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A coupler which integrates surface and subsurface flow implicitly.
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "PK_Helpers.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "mpc_coupled_water.hh"

namespace Amanzi {

const std::string MPCCoupledWater::pk_type_ = "coupled water";

MPCCoupledWater::MPCCoupledWater(const Comm_ptr_type& comm,
                                 Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                                 const Teuchos::RCP<State>& S)
  : MPCStrong<PK_PhysicalBDF_Default>(comm, pk_tree, global_plist, S)
{}


void
MPCCoupledWater::modifyParameterList()
{
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");

  // -- turn on coupling
  pks_list_->sublist(names[0]).set("coupled to surface via flux", true);
  pks_list_->sublist(names[1]).set("coupled to subsurface via flux", true);
  pks_list_->sublist(names[1]).set("scale preconditioner to pressure", true);

  // -- ensure local ops are suface ops
  pks_list_->sublist(names[1]).sublist("diffusion preconditioner").set("surface operator", true);
  pks_list_->sublist(names[1]).sublist("accumulation preconditioner").set("surface operator", true);

  MPCStrong<PK_PhysicalBDF_Default>::modifyParameterList();
}


void
MPCCoupledWater::parseParameterList()
{
  MPCStrong<PK_PhysicalBDF_Default>::parseParameterList();

  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");
  domain_ss_ = pks_list_->sublist(names[0]).get<std::string>("domain name", "domain");
  domain_surf_ = pks_list_->sublist(names[1]).get<std::string>("domain name", "surface");

  // keys
  exfilt_key_ =
    Keys::readKey(*plist_, domain_surf_, "exfiltration flux", "surface_subsurface_flux");
}


void
MPCCoupledWater::setup()
{
  // grab the meshes
  surf_mesh_ = S_->GetMesh(domain_surf_);
  domain_mesh_ = S_->GetMesh(domain_ss_);

  // cast the PKs
  domain_flow_pk_ = sub_pks_[0];
  surf_flow_pk_ = sub_pks_[1];

  // call the MPC's setup, which calls the sub-pk's setups
  MPCStrong<PK_PhysicalBDF_Default>::setup();

  // require the coupling fields, claim ownership
  S_->Require<CompositeVector, CompositeVectorSpace>(exfilt_key_, tag_next_, exfilt_key_)
    .SetMesh(surf_mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  PKHelpers::requireEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // Create the preconditioner.
  // -- collect the preconditioners
  precon_ = domain_flow_pk_->getPreconditioner();
  precon_surf_ = surf_flow_pk_->getPreconditioner();

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
  water_->setTags(tag_current_, tag_next_);
  water_->set_indices(0, 1);

  // grab the debuggers
  domain_db_ = domain_flow_pk_->getDebugger();
  water_->set_db(domain_db_);
  surf_db_ = surf_flow_pk_->getDebugger();

  // With this MPC, thanks to the form of the error/solver, it is often easier
  // to figure out what subsurface face or cell to debug, and it is hard to
  // figure out what the corresponding cell of the surface system is.
  // Therefore, for all debug cells of the subsurface, if that cell is in the
  // top layer of cells, we add the corresponding face's surface cell.
  AmanziMesh::Entity_ID_List debug_cells = domain_db_->get_cells();

  int ncells_surf =
    surf_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  if (debug_cells.size() > 0) {
    const auto& domain_cell_map = domain_mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
    const auto& surf_cell_map = surf_mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);

    std::vector<AmanziMesh::Entity_ID> surf_debug_cells;

    for (int sc = 0; sc != ncells_surf; ++sc) {
      int f = surf_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      auto c = AmanziMesh::getFaceOnBoundaryInternalCell(*domain_mesh_, f);

      auto c_gid = domain_mesh_->getEntityGID(AmanziMesh::Entity_kind::CELL, c);
      if (std::find(debug_cells.begin(), debug_cells.end(), c_gid) != debug_cells.end())
        surf_debug_cells.emplace_back(
          surf_mesh_->getEntityGID(AmanziMesh::Entity_kind::CELL, sc));
    }
    if (surf_debug_cells.size() > 0) surf_db_->add_cells(surf_debug_cells);
  }
}


void
MPCCoupledWater::initialize()
{
  // initialize coupling terms
  S_->GetPtrW<CompositeVector>(exfilt_key_, tag_next_, exfilt_key_)->putScalar(0.);
  S_->GetRecordW(exfilt_key_, tag_next_, exfilt_key_).set_initialized();
  PKHelpers::changedEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // this is manually managed to get the order right
  domain_flow_pk_->initialize();

  // ensure continuity of ICs... subsurface takes precedence.
  MPCHelpers::copySubsurfaceToSurface(
    S_->Get<CompositeVector>(Keys::getKey(domain_ss_, "pressure"), tag_next_),
    S_->GetW<CompositeVector>(
      Keys::getKey(domain_surf_, "pressure"), tag_next_, sub_pks_[1]->getName()));

  // let surface override it -- particularly to get faces right
  surf_flow_pk_->initialize();

  // now we can initialize the bdf time integrator with the initial solution
  PK_BDF_Default::initialize();
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
  moveSolutionToState(*u_new, tag_next_);

  // Evaluate the surface flow residual
  surf_flow_pk_->FunctionalResidual(
    t_old, t_new, u_old->getSubVector(1), u_new->getSubVector(1), g->getSubVector(1));

  // The residual of the surface flow equation provides the water flux from
  // subsurface to surface.
  {
    auto bc =
      S_->GetW<CompositeVector>(exfilt_key_, tag_next_, exfilt_key_).viewComponent("cell", false);
    auto g_surf = g->getSubVector(1)->getData()->viewComponent("cell", false);

    // take off the face area factor to allow it to be used as boundary condition
    const AmanziMesh::MeshCache& mc = domain_mesh_->getCache();
    const AmanziMesh::MeshCache& mc_surf = surf_mesh_->getCache();
    Kokkos::parallel_for(
      "MPCCoupledWater::FunctionalResidual copy to BC",
      g_surf.extent(0),
      KOKKOS_LAMBDA(const int sc) {
        AmanziMesh::Entity_ID f = mc_surf.getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        bc(sc, 0) = g_surf(sc, 0) / mc.getFaceArea(f);
      });
  }
  PKHelpers::changedEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->FunctionalResidual(
    t_old, t_new, u_old->getSubVector(0), u_new->getSubVector(0), g->getSubVector(0));

  // All surface to subsurface fluxes have been taken by the subsurface.
  g->getSubVector(1)->getData()->getComponent("cell", false)->putScalar(0.);
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
  int ierr = precon_->applyInverse(*u->getSubVector(0)->getData(), *Pu->getSubVector(0)->getData());

  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying  CopySubsurfaceToSurface." << std::endl;
  // Copy subsurface face corrections to surface cell corrections
  MPCHelpers::copySubsurfaceToSurface(*Pu->getSubVector(0)->getData(),
                                      *Pu->getSubVector(1)->getData());

  // // Derive surface face corrections.
  // UpdateConsistentFaceCorrectionWater_(u, Pu);

  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");

    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(u->getSubVector(0)->getData().ptr());
    vecs.push_back(Pu->getSubVector(0)->getData().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = u->getSubVector(1)->getData().ptr();
    vecs[1] = Pu->getSubVector(1)->getData().ptr();

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
  modified = domain_flow_pk_->ModifyPredictor(h, u0->getSubVector(0), u->getSubVector(0));

  // Merge surface cells with subsurface faces
  if (modified) {
    // -- not used... something smells... --ETC
    //S_->GetEvaluator(Keys::getKey(domain_surf_,"relative_permeability"), tag_next_)->Update(*S_, name_);
    Teuchos::RCP<const CompositeVector> h_prev =
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_surf_, "ponded_depth"), tag_next_);
    MPCHelpers::mergeSubsurfaceAndSurfacePressure(
      *h_prev, *u->getSubVector(0)->getData(), *u->getSubVector(1)->getData());
  }

  // Hack surface faces
  bool newly_modified = false;
  newly_modified |= water_->ModifyPredictor_Heuristic(h, u);
  newly_modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
  modified |= newly_modified;

  // -- copy surf --> sub
  if (newly_modified) {
    MPCHelpers::copySurfaceToSubsurface(*u->getSubVector(1)->getData(),
                                        *u->getSubVector(0)->getData());
  }

  // Calculate consistent surface faces
  modified |= surf_flow_pk_->ModifyPredictor(h, u0->getSubVector(1), u->getSubVector(1));
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
    vnames.push_back("res");
    vnames.push_back("PC*res");

    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(res->getSubVector(0)->getData().ptr());
    vecs.push_back(du->getSubVector(0)->getData().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = res->getSubVector(1)->getData().ptr();
    vecs[1] = du->getSubVector(1)->getData().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);
  }

  // modify correction using sub-pk approaches
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult modified_res =
    MPCStrong<PK_PhysicalBDF_Default>::ModifyCorrection(h, res, u, du);

  // modify correction using water approaches
  int n_modified = 0;
  n_modified += water_->ModifyCorrection_WaterFaceLimiter(h, res, u, du);

  double damping1 = water_->ModifyCorrection_WaterSpurtDamp(h, res, u, du);
  double damping2 = water_->ModifyCorrection_DesaturatedSpurtDamp(h, res, u, du);
  double damping = std::min(damping1, damping2);
  n_modified += water_->ModifyCorrection_WaterSpurtCap(h, res, u, du, damping);
  n_modified += water_->ModifyCorrection_DesaturatedSpurtCap(h, res, u, du, damping);

  // -- accumulate globally
  int n_modified_l = n_modified;
  Teuchos::reduceAll(*domain_mesh_->getComm(), Teuchos::REDUCE_SUM, 1, &n_modified_l, &n_modified);
  bool modified = (n_modified > 0) || (damping < 1.);

  // -- calculate consistent subsurface cells
  if (modified) {
    // if (consistent_cells_) {
    //   // Derive subsurface cell corrections.
    //   precon_->UpdateConsistentCellCorrection(
    //       *u->getSubVector(0)->getData(),
    //       du->getSubVector(0)->getData().ptr());
    // }

    // Copy subsurface face corrections to surface cell corrections
    MPCHelpers::copySubsurfaceToSurface(*du->getSubVector(0)->getData(),
                                        *du->getSubVector(1)->getData());
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
    vecs.push_back(res->getSubVector(0)->getData().ptr());
    vecs.push_back(du->getSubVector(0)->getData().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = res->getSubVector(1)->getData().ptr();
    vecs[1] = du->getSubVector(1)->getData().ptr();

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
  auto res2 = Teuchos::rcp(new TreeVector(*res, Teuchos::DataAccess::Copy, InitMode::COPY));

  {
    auto res_face = res2->getSubVector(0)->getData()->viewComponent("face", false);
    auto res_surf_cell = res2->getSubVector(1)->getData()->viewComponent("cell", false);
    const auto u_surf_cell = u->getSubVector(1)->getData()->viewComponent("cell", false);
    double p_atm = S_->Get<double>("atmospheric_pressure", Tags::NEXT);
    const AmanziMesh::MeshCache& mc_surf = surf_mesh_->getCache();

    Kokkos::parallel_for(
      "MPCCoupledWater::ErrorNorm", u_surf_cell.extent(0), KOKKOS_LAMBDA(const int& c) {
        if (u_surf_cell(c, 0) > p_atm) {
          auto f = mc_surf.getEntityParent(AmanziMesh::Entity_kind::CELL, c);
          res_surf_cell(c, 0) = res_face(f, 0);
          res_face(f, 0) = 0.;
        }
      });
  }
  return MPCStrong<PK_PhysicalBDF_Default>::ErrorNorm(u, res2);
}


// void
// MPCCoupledWater::UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
//         const Teuchos::RCP<TreeVector>& Pu) {
//   Teuchos::OSTab tab = vo_->getOSTab();

//   Teuchos::RCP<CompositeVector> surf_Pp = Pu->getSubVector(1)->getData();
//   Epetra_MultiVector& surf_Pp_c = *surf_Pp->viewComponent("cell",false);

//   Teuchos::RCP<const CompositeVector> surf_p = u->getSubVector(1)->getData();

//   // Calculate delta h on the surface
//   Teuchos::RCP<CompositeVector> surf_Ph = Teuchos::rcp(new CompositeVector(*surf_Pp));
//   surf_Ph->putScalar(0.);

//   // old ponded depth
//   S_next_->GetEvaluator("ponded_depth")->Update(S_next_.ptr(), name_);
//   *surf_Ph->viewComponent("cell",false) = *S_next_->Get<CompositeVector>("ponded_depth").viewComponent("cell",false);

//   // new ponded depth
//   Teuchos::RCP<TreeVector> tv_p = Teuchos::rcp(new TreeVector());
//   Teuchos::RCP<CompositeVector> cv_p = S_next_->GetPtrW<CompositeVector>("surface-pressure", sub_pks_[1]->getName());
//   cv_p->viewComponent("cell",false)->Update(-1., surf_Pp_c, 1.);
//   tv_p->SetData(cv_p);

//   sub_pks_[1]->markChangedSolution();

//   if (sub_pks_[1]->IsAdmissible(tv_p)) {
//     S_next_->GetEvaluator("ponded_depth")->Update(S_next_.ptr(), name_);

//     // put delta ponded depth into surf_Ph_cell
//     surf_Ph->viewComponent("cell",false)
//         ->Update(-1., *S_next_->Get<CompositeVector>("ponded_depth").viewComponent("cell",false), 1.);

//     // update delta faces
//     precon_surf_->UpdateConsistentFaceCorrection(*surf_p, surf_Ph.ptr());
//     *surf_Pp->viewComponent("face",false) = *surf_Ph->viewComponent("face",false);

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
//   S_next_->GetPtrW<CompositeVector>("surface-pressure",sub_pks_[1]->getName())
//       ->viewComponent("cell",false)->Update(1., surf_Pp_c, 1.);
//   sub_pks_[1]->markChangedSolution();
// }


} // namespace Amanzi
