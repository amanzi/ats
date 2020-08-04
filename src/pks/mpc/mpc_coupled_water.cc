#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "LinearOperatorFactory.hh"

#include "mpc_surface_subsurface_helpers.hh"
#include "mpc_coupled_water.hh"

namespace Amanzi {

MPCCoupledWater::MPCCoupledWater(Teuchos::ParameterList& FElist,
                  const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln) :
    PK(FElist, plist,  S, soln),
    StrongMPC<PK_PhysicalBDF_Default>(FElist, plist,  S, soln) {}



void
MPCCoupledWater::Setup(const Teuchos::Ptr<State>& S) {
  // tweak the sub-PK parameter lists
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");

  // -- turn on coupling
  pks_list_->sublist(names[0]).set("coupled to surface via flux", true);
  pks_list_->sublist(names[1]).set("coupled to subsurface via flux", true);

  // -- ensure local ops are suface ops
  pks_list_->sublist(names[1]).sublist("diffusion preconditioner").set("surface operator", true);
  pks_list_->sublist(names[1]).sublist("accumulation preconditioner").set("surface operator", true);

  
  domain_ss_ = plist_->get<std::string>("domain name","domain");
  domain_surf_ = (domain_ss_.empty() || domain_ss_ == "domain") ? "surface" : std::string("surface_")+domain_ss_;
  domain_surf_ = plist_->get<std::string>("surface domain name",domain_surf_);
  // grab the meshes 
  surf_mesh_ = S->GetMesh(domain_surf_);
  domain_mesh_ = S->GetMesh(domain_ss_);
  
  // cast the PKs
  domain_flow_pk_ = sub_pks_[0];
  surf_flow_pk_ = sub_pks_[1];
  
  // call the MPC's setup, which calls the sub-pk's setups
  StrongMPC<PK_PhysicalBDF_Default>::Setup(S);

  // require the coupling fields, claim ownership
  S->RequireField(Keys::getKey(domain_surf_,"surface_subsurface_flux"), name_)
    ->SetMesh(surf_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // Create the preconditioner.
  // -- collect the preconditioners
 
  precon_ = domain_flow_pk_->preconditioner();
  precon_surf_ = surf_flow_pk_->preconditioner();

  // -- push the surface local ops into the subsurface global operator
  for (Operators::Operator::op_iterator op = precon_surf_->begin();
       op != precon_surf_->end(); ++op) {
    precon_->OpPushBack(*op);
  }

  // -- must re-symbolic assemble subsurf operators, now that they have a surface operator
  precon_->SymbolicAssembleMatrix();
  precon_->InitializePreconditioner(plist_->sublist("preconditioner"));


  // Potentially create a linear solver
  if (plist_->isSublist("linear solver")) {
    Teuchos::ParameterList linsolve_sublist = plist_->sublist("linear solver");
    AmanziSolvers::LinearOperatorFactory<Operators::Operator,CompositeVector,CompositeVectorSpace> fac;
    lin_solver_ = fac.Create(linsolve_sublist, precon_);
  } else {
    lin_solver_ = precon_;
  }
  
  // set up the Water delegate
  Teuchos::RCP<Teuchos::ParameterList> water_list = Teuchos::sublist(plist_, "water delegate");
  water_ = Teuchos::rcp(new MPCDelegateWater(water_list, domain_ss_));
  water_->set_indices(0,1);

  // grab the debuggers
  domain_db_ = domain_flow_pk_->debugger();
  surf_db_ = surf_flow_pk_->debugger();
  water_->set_db(domain_db_);
}

void

MPCCoupledWater::Initialize(const Teuchos::Ptr<State>& S) {
   // initialize coupling terms
  S->GetFieldData(Keys::getKey(domain_surf_,"surface_subsurface_flux"), name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_surf_,"surface_subsurface_flux"), name_)->set_initialized();
  // Initialize all sub PKs.

  MPC<PK_PhysicalBDF_Default>::Initialize(S);

  // // ensure continuity of ICs... surface takes precedence.
  // CopySurfaceToSubsurface(*S->GetFieldData(Keys::getKey(domain_surf_,"pressure"), sub_pks_[1]->name()),
  //       		  S->GetFieldData(Keys::getKey(domain_ss_,"pressure"), sub_pks_[0]->name()).ptr());
  // ensure continuity of ICs... subsurface takes precedence.
  CopySubsurfaceToSurface(*S->GetFieldData(Keys::getKey(domain_ss_,"pressure"), sub_pks_[0]->name()),
			  S->GetFieldData(Keys::getKey(domain_surf_,"pressure"), sub_pks_[1]->name()).ptr());

  // Initialize my timestepper.
  PK_BDF_Default::Initialize(S);
}

void
MPCCoupledWater::set_states(const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<State>& S_inter,
                            const Teuchos::RCP<State>& S_next) {
  StrongMPC<PK_PhysicalBDF_Default>::set_states(S,S_inter,S_next);
  water_->set_states(S,S_inter,S_next);
}


// -- computes the non-linear functional g = g(t,u,udot)
//    By default this just calls each sub pk FunctionalResidual().
void
MPCCoupledWater::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                            Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // propagate updated info into state
  Solution_to_State(*u_new, S_next_);

  // Evaluate the surface flow residual
  surf_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(1),
                            u_new->SubVector(1), g->SubVector(1));

  // The residual of the surface flow equation provides the mass flux from
  // subsurface to surface.

  Epetra_MultiVector& source = *S_next_->GetFieldData(Keys::getKey(domain_surf_,"surface_subsurface_flux"),
          name_)->ViewComponent("cell",false);
  source = *g->SubVector(1)->Data()->ViewComponent("cell",false);

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0),
          u_new->SubVector(0), g->SubVector(0));

  // All surface to subsurface fluxes have been taken by the subsurface.
  g->SubVector(1)->Data()->ViewComponent("cell",false)->PutScalar(0.);
}

// -- Apply preconditioner to u and returns the result in Pu.
int MPCCoupledWater::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon application:" << std::endl;

  // call the precon's inverse
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying subsurface operator." << std::endl;
  int ierr = lin_solver_->ApplyInverse(*u->SubVector(0)->Data(), *Pu->SubVector(0)->Data());

  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying  CopySubsurfaceToSurface." << std::endl;
  // Copy subsurface face corrections to surface cell corrections
  CopySubsurfaceToSurface(*Pu->SubVector(0)->Data(),
                          Pu->SubVector(1)->Data().ptr());

  // // Derive surface face corrections.
  // UpdateConsistentFaceCorrectionWater_(u, Pu);

  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
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

// -- Update the preconditioner.
void
MPCCoupledWater::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // order important -- subsurface's pk includes the surface's local ops, so
  // doing the subsurface 2nd re-inits the surface matrices (and doesn't
  // refill them).  This is why subsurface is first
  StrongMPC<PK_PhysicalBDF_Default>::UpdatePreconditioner(t, up, h);
  
  precon_->AssembleMatrix();
  precon_->UpdatePreconditioner();
  
}

// -- Modify the predictor.
bool
MPCCoupledWater::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  bool modified = false;

  // Calculate consistent faces
  modified = domain_flow_pk_->ModifyPredictor(h, u0->SubVector(0), u->SubVector(0));

  // Merge surface cells with subsurface faces
  if (modified) {

    S_next_->GetFieldEvaluator(Keys::getKey(domain_surf_,"relative_permeability"))->HasFieldChanged(S_next_.ptr(),name_);
    Teuchos::RCP<const CompositeVector> h_prev = S_inter_->GetFieldData(Keys::getKey(domain_surf_,"ponded_depth"));
    MergeSubsurfaceAndSurfacePressure(*h_prev, u->SubVector(0)->Data().ptr(),
            u->SubVector(1)->Data().ptr());
  }
  

  // Hack surface faces
  bool newly_modified = false;
  newly_modified |= water_->ModifyPredictor_Heuristic(h, u);
  newly_modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
  modified |= newly_modified;

  // -- copy surf --> sub
  if (newly_modified) {
    CopySurfaceToSubsurface(*u->SubVector(1)->Data(), u->SubVector(0)->Data().ptr());
  }

  // Calculate consistent surface faces
  modified |= surf_flow_pk_->ModifyPredictor(h, u0->SubVector(1), u->SubVector(1));

  return modified;
}

// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCCoupledWater::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();
  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "NKA'd, PC'd correction." << std::endl;

    std::vector<std::string> vnames;
    vnames.push_back("p");
    vnames.push_back("PC*p");

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
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
  n_modified += water_->ModifyCorrection_WaterFaceLimiter(h, res, u, du);
  double damping = water_->ModifyCorrection_WaterSpurtDamp(h, res, u, du);
  n_modified += water_->ModifyCorrection_WaterSpurtCap(h, res, u, du, damping);

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
    CopySubsurfaceToSurface(*du->SubVector(0)->Data(),
                            du->SubVector(1)->Data().ptr());
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

    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(res->SubVector(0)->Data().ptr());
    vecs.push_back(du->SubVector(0)->Data().ptr());

    *vo_->os() << " Subsurface precon:" << std::endl;
    domain_db_->WriteVectors(vnames, vecs, true);

    vecs[0] = res->SubVector(1)->Data().ptr();
    vecs[1] = du->SubVector(1)->Data().ptr();

    *vo_->os() << " Surface precon:" << std::endl;
    surf_db_->WriteVectors(vnames, vecs, true);
  }

  return (modified_res || modified) ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING :
      AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
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
//   S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
//   *surf_Ph->ViewComponent("cell",false) = *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false);

//   // new ponded depth
//   Teuchos::RCP<TreeVector> tv_p = Teuchos::rcp(new TreeVector());
//   Teuchos::RCP<CompositeVector> cv_p = S_next_->GetFieldData("surface-pressure", sub_pks_[1]->name());
//   cv_p->ViewComponent("cell",false)->Update(-1., surf_Pp_c, 1.);
//   tv_p->SetData(cv_p);

//   sub_pks_[1]->ChangedSolution();

//   if (sub_pks_[1]->IsAdmissible(tv_p)) {
//     S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);

//     // put delta ponded depth into surf_Ph_cell
//     surf_Ph->ViewComponent("cell",false)
//         ->Update(-1., *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false), 1.);

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
//       vecs.push_back(S_next_->GetFieldData("ponded_depth").ptr());
//       vecs.push_back(surf_Ph.ptr());
//       surf_db_->WriteVectors(vnames, vecs, true);
//     }
//   }

//   // revert solution so we don't break things
//   S_next_->GetFieldData("surface-pressure",sub_pks_[1]->name())
//       ->ViewComponent("cell",false)->Update(1., surf_Pp_c, 1.);
//   sub_pks_[1]->ChangedSolution();
// }


} // namespace
