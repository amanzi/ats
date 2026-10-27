/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "PK_Helpers.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "richards.hh"
#include "overland_pressure.hh"
#include "mpc_coupled_water_noninjected.hh"

namespace Amanzi {

// needs a better name!
const std::string MPCCoupledWaterNonInjected::pk_type_ = "coupled water not injected";

MPCCoupledWaterNonInjected::MPCCoupledWaterNonInjected(
  const Comm_ptr_type& comm,
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
  const Teuchos::RCP<State>& S)
  : MPCStrong<PK_PhysicalBDF_Default>(comm, pk_tree, global_plist, S),
    symbolic_assembled_(false)
{}


void
MPCCoupledWaterNonInjected::parseParameterList()
{
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");
  domain_ss_ = pks_list_->sublist(names[0]).get<std::string>("domain name", "domain");
  domain_surf_ = pks_list_->sublist(names[1]).get<std::string>("domain name", "surface");

  // exfiltration -- flux from subsurface to surface
  exfilt_key_ =
    Keys::readKey(*plist_, domain_surf_, "exfiltration flux", "surface_subsurface_flux");

  // set the coupling flags
  pks_list_->sublist(names[0]).set("coupled to surface via head", true);

  // this flag is named a little weird, but the surface equation gets the
  // Dirichlet-to-Neumann mapped flux from the subsurface
  pks_list_->sublist(names[1]).set("coupled to subsurface via head", true);
  pks_list_->sublist(names[1]).set("surface-subsurface flux key", exfilt_key_);
  pks_list_->sublist(names[1]).set("scale preconditioner to pressure", true);

  // Single switch for all Newton-correction dkr/dp terms (Richards' own
  // diagonal term plus this MPC's coupling-block cross terms) -- read from
  // our own plist and forwarded into the subsurface sub-PK's "diffusion
  // preconditioner" sublist so the two can never disagree.
  include_newton_correction_ = plist_->get<bool>("include Newton correction", false);
  pks_list_->sublist(names[0])
    .sublist("diffusion preconditioner")
    .set("include Newton correction", include_newton_correction_);

  // also need a few keys from other lists
  ss_uw_coef_key_ = Keys::readKey(
    *plist_, domain_ss_, "upwinded conductivity", "upwind_relative_hydraulic_conductivity");
  ss_flux_key_ = Keys::readKey(*plist_, domain_ss_, "darcy flux", "water_flux");
  ss_pressure_key_ = Keys::getKey(domain_ss_, "pressure");

  MPCStrong<PK_PhysicalBDF_Default>::parseParameterList();
}


void
MPCCoupledWaterNonInjected::setup()
{
  // grab the meshes
  surf_mesh_ = S_->GetMesh(domain_surf_);
  domain_mesh_ = S_->GetMesh(domain_ss_);

  // cast the PKs
  domain_flow_pk_ = sub_pks_[0];
  surf_flow_pk_ = sub_pks_[1];

  // require fields
  PKHelpers::requireAtNext(exfilt_key_, tag_next_, *S_, exfilt_key_)
    .SetMesh(surf_mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // the coupling blocks need d(uw_coef)/d(p_sub) on both DoFs (cell and
  // boundary) -- require it only when Newton correction is enabled, same
  // switch Richards itself uses for its own diagonal term.
  if (include_newton_correction_) {
    S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        ss_uw_coef_key_, tag_next_, ss_pressure_key_, tag_next_)
      .SetGhosted();
  }

  // call the MPC's setup, which calls the sub-pk's setups
  MPCStrong<PK_PhysicalBDF_Default>::setup();

  // Build the combined preconditioner as a 2x2 TreeOperator: the two
  // sub-PKs' own preconditioners on the diagonal, with rectangular
  // off-diagonal coupling blocks.
  //
  // -- Build the diagonal block
  auto pc_sub = domain_flow_pk_->getPreconditioner();
  auto pc_surf = surf_flow_pk_->getPreconditioner();
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(pc_sub->getDomainMap())));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(pc_surf->getDomainMap())));
  precon_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  precon_->set_operator_block(0, 0, pc_sub);
  precon_->set_operator_block(1, 1, pc_surf);

  // -- Build the injection map: from surface-cell-local-index sc to
  //    subsurface-cell-global index of its parent subsurface top cell
  int ncells_surf_owned =
    surf_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  Kokkos::View<GO*, Kokkos::HostSpace> parent_gids_host("coupling injection GIDs",
                                                        ncells_surf_owned);
  for (int sc = 0; sc != ncells_surf_owned; ++sc) {
    AmanziMesh::Entity_ID f = surf_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
    AmanziMesh::Entity_ID c = AmanziMesh::getFaceOnBoundaryInternalCell(*domain_mesh_, f);
    parent_gids_host(sc) = domain_mesh_->getEntityGID(AmanziMesh::Entity_kind::CELL, c);
  }
  auto injection = Teuchos::rcp(new Map_type(-1, parent_gids_host, 0, surf_mesh_->getComm()));

  // Our coupling blocks are purely cell-to-cell, regardless of what
  // schema the diagonal blocks (pc_sub/pc_surf) actually use
  Operators::Schema cell_schema(AmanziMesh::Entity_kind::CELL, 1);

  // -- (0,1): row = subsurface cells, col = surface cells
  Teuchos::ParameterList coupling_plist_01("coupling (0,1)");
  coupling_op_01_ = Teuchos::rcp(new Operators::Operator_Schema(
    pc_sub->getRangeMap(), pc_surf->getDomainMap(), coupling_plist_01,
    cell_schema, cell_schema));
  coupling_local_op_01_ = Teuchos::rcp(new Operators::Op_MeshInjection(
    cell_schema, domain_mesh_, cell_schema, surf_mesh_, injection,
    /*transpose=*/true));
  coupling_op_01_->OpPushBack(coupling_local_op_01_);
  precon_->set_operator_block(0, 1, coupling_op_01_);

  // -- (1,0): row = surface cells, col = subsurface cells
  Teuchos::ParameterList coupling_plist_10("coupling (1,0)");
  coupling_op_10_ = Teuchos::rcp(new Operators::Operator_Schema(
    pc_surf->getRangeMap(), pc_sub->getDomainMap(), coupling_plist_10,
    cell_schema, cell_schema));
  coupling_local_op_10_ = Teuchos::rcp(new Operators::Op_MeshInjection(
    cell_schema, surf_mesh_, cell_schema, domain_mesh_, injection,
    /*transpose=*/false));
  coupling_op_10_->OpPushBack(coupling_local_op_10_);
  precon_->set_operator_block(1, 0, coupling_op_10_);

  // -- (1,1) diagonal surface, d(exfilt)/dp_surf
  coupling_acc_surf_ = Teuchos::rcp(
    new Operators::Op_Cell_Cell("coupling (1,1) diagonal addition", surf_mesh_));
  pc_surf->OpPushBack(coupling_acc_surf_);

  // -- set parameters for an inverse
  precon_->set_inverse_parameters(plist_->sublist("inverse"));

  // set up the Water delegate
  Teuchos::RCP<Teuchos::ParameterList> water_list = Teuchos::sublist(plist_, "water delegate");
  water_ = Teuchos::rcp(new MPCDelegateWater(water_list, S_, domain_ss_, domain_surf_));
  water_->setTags(tag_current_, tag_next_);
  water_->set_indices(0, 1);

  // grab the debuggers
  domain_db_ = domain_flow_pk_->getDebugger();
  water_->set_db(domain_db_);
  surf_db_ = surf_flow_pk_->getDebugger();
}


void
MPCCoupledWaterNonInjected::initialize()
{
  // initialize the exchange-flux diagnostic
  S_->GetPtrW<CompositeVector>(exfilt_key_, tag_next_, exfilt_key_)->putScalar(0.);
  S_->GetRecordW(exfilt_key_, tag_next_, exfilt_key_).set_initialized();
  PKHelpers::changedEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // this is manually managed to get the order right
  domain_flow_pk_->initialize();

  // ensure continuity of ICs
  MPCHelpers::copySubsurfaceToSurface(
    S_->Get<CompositeVector>(Keys::getKey(domain_ss_, "pressure"), tag_next_),
    S_->GetW<CompositeVector>(
      Keys::getKey(domain_surf_, "pressure"), tag_next_, sub_pks_[1]->getName()));

  // let surface override it -- particularly to get faces right
  surf_flow_pk_->initialize();

  // now we can initialize the bdf time integrator with the initial solution
  PK_BDF_Default::initialize();
}


// -- computes the non-linear functional g = g(t,u,udot).
void
MPCCoupledWaterNonInjected::FunctionalResidual(double t_old,
                                               double t_new,
                                               Teuchos::RCP<TreeVector> u_old,
                                               Teuchos::RCP<TreeVector> u_new,
                                               Teuchos::RCP<TreeVector> g)
{
  moveSolutionToState(*u_new, tag_next_);

  // 1. Subsurface residual, with surface pressure imposed as a Dirichlet BC
  domain_flow_pk_->FunctionalResidual(
    t_old, t_new, u_old->getSubVector(0), u_new->getSubVector(0), g->getSubVector(0));

  // 2/3. The subsurface's ApplyDiffusion_/UpdateFlux populated the flux field
  //      on all faces, including the top boundary face, consistent with the
  //      Dirichlet BC used in step 1 -- stuff that face's flux into
  //      exfiltration, oriented outward from the subsurface domain
  {
    auto exfilt =
      S_->GetW<CompositeVector>(exfilt_key_, tag_next_, exfilt_key_).viewComponent("cell", false);
    const auto flux =
      S_->Get<CompositeVector>(ss_flux_key_, tag_next_).viewComponent("face", false);
    const AmanziMesh::MeshCache& mc = domain_mesh_->getCache();
    const AmanziMesh::MeshCache& mc_surf = surf_mesh_->getCache();
    Kokkos::parallel_for(
      "MPCCoupledWaterNonInjected::FunctionalResidual exfiltration flux",
      exfilt.extent(0),
      KOKKOS_LAMBDA(const int sc) {
        AmanziMesh::Entity_ID f = mc_surf.getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        auto cells = mc.getFaceCells(f);
        int c = cells(0);
        auto [cfaces, dirs] = mc.getCellFacesAndDirections(c);
        double dir = 1.0;
        for (int n = 0; n != cfaces.extent(0); ++n) {
          if (cfaces(n) == f) {
            dir = dirs(n);
            break;
          }
        }
        // flux(f) is signed outward-positive relative to cell c when dir>0
        exfilt(sc, 0) = dir * flux(f, 0);
      });
  }
  PKHelpers::changedEvaluatorPrimary(exfilt_key_, tag_next_, *S_);

  // 4. Surface residual, consuming exfilt_key_ as its
  //    coupled-to-subsurface-via-head source term.
  surf_flow_pk_->FunctionalResidual(
    t_old, t_new, u_old->getSubVector(1), u_new->getSubVector(1), g->getSubVector(1));
}


void
MPCCoupledWaterNonInjected::UpdatePreconditioner(double t,
                                                 Teuchos::RCP<const TreeVector> up,
                                                 double h)
{
  // base (0,0) and (1,1) blocks computed from base PKs
  MPCStrong<PK_PhysicalBDF_Default>::UpdatePreconditioner(t, up, h);

  // (0,1)/(1,0) off-diagonal blocks, plus the (1,1)-diagonal addition for
  // overland flow d(exfilt)/dp_surf
  UpdateCouplingBlocks_();

  // manually assembled until the bug fix in Amanzi master is propagated to
  // this branch
  if (!symbolic_assembled_) {
    precon_->SymbolicAssembleMatrix();
    symbolic_assembled_ = true;
  }
  precon_->AssembleMatrix();
}


void
MPCCoupledWaterNonInjected::UpdateCouplingBlocks_()
{
  Teuchos::RCP<Flow::Richards> richards =
    Teuchos::rcp_dynamic_cast<Flow::Richards>(domain_flow_pk_);
  AMANZI_ASSERT(richards != Teuchos::null);

  auto getTransGrav = richards->getDiffusionPreconditioner()->GetTransmissibilityGravityFunctor();

  const auto uw_coef_f =
    S_->Get<CompositeVector>(ss_uw_coef_key_, tag_next_).viewComponent("face", false);

  // Newton-correction dkr/dp terms -- all on or all off, same switch
  // Richards uses for its own diagonal term (see parseParameterList()).
  // When off, duw_coef_f is never fetched (its derivative was never
  // Require()d either) and the kernel below treats both DoFs as zero.
  bool include_newton_correction = include_newton_correction_;
  const auto duw_coef_f = include_newton_correction_ ?
    S_->GetDerivative<CompositeVector>(ss_uw_coef_key_, tag_next_, ss_pressure_key_, tag_next_)
      .viewComponent("face", false) :
    decltype(S_->GetDerivative<CompositeVector>(
                 ss_uw_coef_key_, tag_next_, ss_pressure_key_, tag_next_)
               .viewComponent("face", false))();

  const auto p_sub_c = S_->Get<CompositeVector>(ss_pressure_key_, tag_next_).viewComponent("cell", false);
  const auto p_surf_c =
    S_->Get<CompositeVector>(Keys::getKey(domain_surf_, "pressure"), tag_next_)
      .viewComponent("cell", false);

  const AmanziMesh::MeshCache& mc_surf = surf_mesh_->getCache();
  const AmanziMesh::MeshCache& mc = domain_mesh_->getCache();

  auto diag01_v = coupling_local_op_01_->diag->getLocalViewDevice(Tpetra::Access::ReadWrite);
  auto diag10_v = coupling_local_op_10_->diag->getLocalViewDevice(Tpetra::Access::ReadWrite);
  auto diag11_v = coupling_acc_surf_->diag->getLocalViewDevice(Tpetra::Access::ReadWrite);

  int ncells_surf = diag01_v.extent(0);

  Kokkos::parallel_for(
    "MPCCoupledWaterNonInjected::UpdateCouplingBlocks_",
    ncells_surf,
    KOKKOS_LAMBDA(const int sc) {
      AmanziMesh::Entity_ID f = mc_surf.getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      auto cells = mc.getFaceCells(f);
      int c = cells(0);

      // outward orientation of f relative to the subsurface cell c, same
      // idiom used for the exfiltration flux sign in FunctionalResidual.
      auto [cfaces, dirs] = mc.getCellFacesAndDirections(c);
      double dir = 1.0;
      for (int n = 0; n != cfaces.extent(0); ++n) {
        if (cfaces(n) == f) {
          dir = dirs(n);
          break;
        }
      }

      auto trans_grav_f = getTransGrav(f);
      double trans = trans_grav_f.first;
      double grav = trans_grav_f.second;
      double uw_coef = uw_coef_f(f, 0);
      double dpres = p_sub_c(c, 0) - p_surf_c(sc, 0);

      // q(p_sub, p_surf) = trans * uw_coef * dpres + dir * uw_coef * grav,
      // the Dirichlet-face flux (same formula as PDE_DiffusionFVwithGravity
      // ::AnalyticJacobian_'s own boundary-face term). duw_coef_f(f,0/1)
      // are d(uw_coef)/d(p_sub) and d(uw_coef)/d(p_surf) respectively --
      // left at 0 (the frozen-mobility approximation) when Newton
      // correction is disabled.
      double duw_dp_sub = include_newton_correction ? duw_coef_f(f, 0) : 0.0;
      double duw_dp_surf = include_newton_correction ? duw_coef_f(f, 1) : 0.0;

      double dq_dp_sub = trans * uw_coef + (trans * dpres + dir * grav) * duw_dp_sub;
      double dq_dp_surf = -trans * uw_coef + (trans * dpres + dir * grav) * duw_dp_surf;

      // (0,1): d(g_sub)/d(p_surf) == d(q)/d(p_surf).
      diag01_v(sc, 0) = dq_dp_surf;

      // (1,0): d(g_surf)/d(p_sub) == -d(q)/d(p_sub).
      diag10_v(sc, 0) = -dq_dp_sub;

      // (1,1) diagonal ADDITION: d(g_surf)/d(p_surf) == -d(q)/d(p_surf).
      diag11_v(sc, 0) += -dq_dp_surf;
    });
}


int
MPCCoupledWaterNonInjected::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                                Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Precon application:" << std::endl;

  int ierr = precon_->applyInverse(*u, *Pu);

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
MPCCoupledWaterNonInjected::ModifyPredictor(double h,
                                            Teuchos::RCP<const TreeVector> u0,
                                            Teuchos::RCP<TreeVector> u)
{
  // potentially modify subsurface, if this changes boundary faces, copy over to surface
  bool modified = domain_flow_pk_->ModifyPredictor(h, u0->getSubVector(0), u->getSubVector(0));
  if (modified) {
    MPCHelpers::copySubsurfaceToSurface(*u->getSubVector(0)->getData(),
                                        *u->getSubVector(1)->getData());
  }

  // Water delegate heuristics: read/write the subsurface boundary_face
  // value directly, and (per their own implementation) already write the
  // corresponding surface cell value themselves as their last step -- see
  // MPCDelegateWater::ModifyPredictor_Heuristic/_WaterSpurtDamp -- so no
  // additional sync is needed after this block.
  bool newly_modified = false;
  newly_modified |= water_->ModifyPredictor_Heuristic(h, u);
  newly_modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
  modified |= newly_modified;

  // Surface's own predictor modification changes the surface cell value, and
  // sync with subsurface
  bool surf_modified = surf_flow_pk_->ModifyPredictor(h, u0->getSubVector(1), u->getSubVector(1));
  modified |= surf_modified;
  if (surf_modified) {
    MPCHelpers::copySurfaceToSubsurface(*u->getSubVector(1)->getData(),
                                        *u->getSubVector(0)->getData());
  }

  return modified;
}


// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCCoupledWaterNonInjected::ModifyCorrection(double h,
                                             Teuchos::RCP<const TreeVector> res,
                                             Teuchos::RCP<const TreeVector> u,
                                             Teuchos::RCP<TreeVector> du)
{
  // modify subsurfafce and sync
  auto ss_modified_res = domain_flow_pk_->ModifyCorrection(h, res->getSubVector(0),
          u->getSubVector(0), du->getSubVector(0));
  if (ss_modified_res) {
    MPCHelpers::copySubsurfaceToSurface(*du->getSubVector(0)->getData(),
            *du->getSubVector(1)->getData());
  }

  // modify surface and sync
  auto surf_modified_res = surf_flow_pk_->ModifyCorrection(h, res->getSubVector(1),
          u->getSubVector(1), du->getSubVector(1));
  if (surf_modified_res) {
    MPCHelpers::copySurfaceToSubsurface(*du->getSubVector(1)->getData(),
            *du->getSubVector(0)->getData());
  }

  // Water delegate corrections: read/write the subsurface boundary_face
  // correction only.
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

  // Sync the (possibly capped/damped) subsurface boundary_face correction
  // back into the surface cell correction, so the surface equation's
  // actual solved correction reflects any limiting the delegate applied.
  if (modified) {
    MPCHelpers::copySubsurfaceToSurface(*du->getSubVector(0)->getData(),
                                        *du->getSubVector(1)->getData());
  }

  return (ss_modified_res || surf_modified_res || modified) ?
           AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING :
           AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

} // namespace Amanzi
