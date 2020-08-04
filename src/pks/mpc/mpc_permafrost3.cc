#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "mpc_delegate_ewc_surface.hh"
#include "mpc_delegate_ewc_subsurface.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "permafrost_model.hh"
#include "surface_ice_model.hh"
#include "energy_base.hh"
#include "advection.hh"
#include "LinearOperatorFactory.hh"
#include "mpc_permafrost3.hh"

namespace Amanzi {

MPCPermafrost3::MPCPermafrost3(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
  StrongMPC<PKPhysicalBDFBase>(plist, FElist, soln) {}


void
MPCPermafrost3::setup(const Teuchos::Ptr<State>& S) {
  // tweak the sub-PK parameter lists
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");

  // -- turn off PC assembly in the subsurface unless we need the forward operator
  if (plist_->isSublist("Coupled Solver")) {
    pks_list_->sublist(names[0]).set("assemble preconditioner", true);
    pks_list_->sublist(names[1]).set("assemble preconditioner", true);
  } else {
    pks_list_->sublist(names[0]).set("assemble preconditioner", false);
    pks_list_->sublist(names[1]).set("assemble preconditioner", false);
  }

  // -- always use the TPFA forward operator for surface
  pks_list_->sublist(names[2]).set("assemble preconditioner", true);
  pks_list_->sublist(names[3]).set("assemble preconditioner", true);

  // -- turn on coupling
  pks_list_->sublist(names[0]).set("coupled to surface via flux", true);
  pks_list_->sublist(names[1]).set("coupled to surface via flux", true);
  pks_list_->sublist(names[2]).set("coupled to subsurface via flux", true);
  pks_list_->sublist(names[3]).set("coupled to subsurface via flux", true);

  // -- set up PC coupling
  //  pks_list_->sublist(names[0]).sublist("diffusion preconditioner").set("coupled to surface", true);
  //  pks_list_->sublist(names[1]).sublist("diffusion preconditioner").set("coupled to surface", true);
  pks_list_->sublist(names[2]).sublist("diffusion preconditioner").set("TPFA", true);
  pks_list_->sublist(names[3]).sublist("diffusion preconditioner").set("TPFA", true);

  // grab the meshes
  surf_mesh_ = S->GetMesh("surface");
  domain_mesh_ = S->GetMesh();

  // cast the PKs
  domain_flow_pk_ = sub_pks_[0];
  domain_energy_pk_ = sub_pks_[1];
  surf_flow_pk_ = sub_pks_[2];
  surf_energy_pk_ = sub_pks_[3];

  // call the MPC's setup, which calls the sub-pk's setups
  StrongMPC<PKPhysicalBDFBase>::setup(S);

  // require the coupling fields, claim ownership
  S->RequireField("surface_subsurface_flux", name_)
      ->SetMesh(surf_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("surface_subsurface_energy_flux", name_)
      ->SetMesh(surf_mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // Create the preconditioner.
  // -- Allocate matrix.
  Teuchos::ParameterList& pc_sublist = plist_->sublist("Coupled PC");
  precon_ = Teuchos::rcp(new Operators::MatrixMFD_Coupled_Surf(pc_sublist, domain_mesh_));

  // -- Collect the sub-blocks.
  /*
    pc_flow_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_Surf>(
    domain_flow_pk_->preconditioner());
    ASSERT(pc_flow_ != Teuchos::null);
    pc_energy_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_Surf>(
    domain_energy_pk_->preconditioner());
    ASSERT(pc_energy_ != Teuchos::null);
  */

  pc_flow_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(
      domain_flow_pk_->preconditioner());
  ASSERT(pc_flow_ != Teuchos::null);
  pc_energy_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(
      domain_energy_pk_->preconditioner());
  ASSERT(pc_energy_ != Teuchos::null);

  pc_surf_flow_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(
      surf_flow_pk_->preconditioner());
  ASSERT(pc_surf_flow_ != Teuchos::null);
  pc_surf_energy_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(
      surf_energy_pk_->preconditioner());
  ASSERT(pc_surf_energy_ != Teuchos::null);

  /*
  // -- Subsurface blocks include their surface operators.
  pc_flow_->SetSurfaceOperator(pc_surf_flow_);
  pc_energy_->SetSurfaceOperator(pc_surf_energy_);

  // -- must re-symbolic assemble surf operators, now that they have a surface operator
  pc_flow_->SymbolicAssembleGlobalMatrices();
  pc_energy_->SymbolicAssembleGlobalMatrices();
  */

  // -- finally, set the sub blocks in the coupled PC
  precon_->SetSubBlocks(pc_flow_, pc_energy_);
  precon_->SetSurfaceOperators(pc_surf_flow_, pc_surf_energy_);
  precon_->SymbolicAssembleGlobalMatrices();
  precon_->InitPreconditioner();

  // set up the advective term
  // -- clone the flow operator
  Teuchos::RCP<CompositeMatrix> pcAdv_mat = sub_pks_[0]->preconditioner()->Clone();
  pcAdv_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(pcAdv_mat);
  ASSERT(pcAdv_ != Teuchos::null);
  precon_->SetAdvectiveBlock(pcAdv_);
  // -- get the field -- this is very hackish and demonstrates why coupled PKs should be redesigned --etc
  Teuchos::RCP<Energy::EnergyBase> pk_as_energy = Teuchos::rcp_dynamic_cast<Energy::EnergyBase>(domain_energy_pk_);
  ASSERT(pk_as_energy != Teuchos::null);
  adv_field_ = pk_as_energy->advection()->field();
  adv_flux_ = pk_as_energy->advection()->flux();
  
  // Potential create the solver
  if (plist_->isSublist("Coupled Solver")) {
    Teuchos::ParameterList linsolve_sublist = plist_->sublist("Coupled Solver");
    AmanziSolvers::LinearOperatorFactory<TreeMatrix,TreeVector,TreeVectorSpace> fac;
    lin_solver_ = fac.Create(linsolve_sublist, precon_);
  } else {
    lin_solver_ = precon_;
  }

  // select the method used for preconditioning
  std::string precon_string = plist_->get<std::string>("preconditioner type", "picard");
  if (precon_string == "none") {
    precon_type_ = PRECON_NONE;
  } else if (precon_string == "picard") {
    precon_type_ = PRECON_PICARD;
  } else if (precon_string == "ewc") {
    precon_type_ = PRECON_EWC;
  } else if (precon_string == "smart ewc") {
    precon_type_ = PRECON_EWC;
  } else {
    Errors::Message message(std::string("Invalid preconditioner type ")+precon_string);
    Exceptions::amanzi_throw(message);
  }

  // set up the EWC delegates
  Teuchos::RCP<Teuchos::ParameterList> sub_ewc_list = Teuchos::sublist(plist_, "subsurface ewc delegate");
  sub_ewc_list->set("PK name", name_);
  sub_ewc_list->set("domain key", "");
  sub_ewc_ = Teuchos::rcp(new MPCDelegateEWCSubsurface(*sub_ewc_list));
  Teuchos::RCP<PermafrostModel> sub_model = Teuchos::rcp(new PermafrostModel());
  sub_ewc_->set_model(sub_model);
  sub_ewc_->setup(S);

  Teuchos::RCP<Teuchos::ParameterList> surf_ewc_list = Teuchos::sublist(plist_, "surface ewc delegate");
  surf_ewc_list->set("PK name", name_);
  surf_ewc_list->set("domain key", "surface");
  surf_ewc_ = Teuchos::rcp(new MPCDelegateEWCSurface(*surf_ewc_list));
  Teuchos::RCP<SurfaceIceModel> surf_model = Teuchos::rcp(new SurfaceIceModel());
  surf_ewc_->set_model(surf_model);
  surf_ewc_->setup(S);

  // set up the Water delegate
  Teuchos::RCP<Teuchos::ParameterList> water_list = Teuchos::sublist(plist_, "water delegate");
  water_ = Teuchos::rcp(new MPCDelegateWater(water_list));
  water_->set_indices(0,2,1,3);

  // set up our own predictors
  consistent_cells_ =
      plist_->get<bool>("ensure consistent cells after face updates", false);

  // grab the debuggers
  domain_db_ = domain_flow_pk_->debugger();
  surf_db_ = surf_flow_pk_->debugger();
}

void
MPCPermafrost3::initialize(const Teuchos::Ptr<State>& S) {
  // initialize coupling terms
  S->GetFieldData("surface_subsurface_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_flux", name_)->set_initialized();
  S->GetFieldData("surface_subsurface_energy_flux", name_)->PutScalar(0.);
  S->GetField("surface_subsurface_energy_flux", name_)->set_initialized();

  // Initialize all sub PKs.
  MPC<PKPhysicalBDFBase>::initialize(S);

  // ensure continuity of ICs... surface takes precedence.
  CopySurfaceToSubsurface(*S->GetFieldData("surface_pressure", sub_pks_[2]->name()),
                          S->GetFieldData("pressure", sub_pks_[0]->name()).ptr());
  CopySurfaceToSubsurface(*S->GetFieldData("surface_temperature", sub_pks_[3]->name()),
                          S->GetFieldData("temperature", sub_pks_[1]->name()).ptr());

  // initialize delegates
  surf_ewc_->initialize(S);
  sub_ewc_->initialize(S);

  // Initialize my timestepper.
  PKBDFBase::initialize(S);

  // // advection mass matrices
  // *pcAdv_ = *pc_flow_;
}

void
MPCPermafrost3::set_states(const Teuchos::RCP<const State>& S,
                           const Teuchos::RCP<State>& S_inter,
                           const Teuchos::RCP<State>& S_next) {
  StrongMPC<PKPhysicalBDFBase>::set_states(S,S_inter,S_next);
  surf_ewc_->set_states(S,S_inter,S_next);
  sub_ewc_->set_states(S,S_inter,S_next);
  water_->set_states(S,S_inter,S_next);
}


void
MPCPermafrost3::commit_state(double dt, const Teuchos::RCP<State>& S) {
  StrongMPC<PKPhysicalBDFBase>::commit_state(dt,S);
  if (surf_ewc_ != Teuchos::null) surf_ewc_->commit_state(dt,S);
  if (sub_ewc_ != Teuchos::null) sub_ewc_->commit_state(dt,S);
}


// -- computes the non-linear functional g = g(t,u,udot)
//    By default this just calls each sub pk Functional().
void
MPCPermafrost3::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // propagate updated info into state
  solution_to_state(*u_new, S_next_);

  // Evaluate the surface flow residual
  surf_flow_pk_->Functional(t_old, t_new, u_old->SubVector(2),
                            u_new->SubVector(2), g->SubVector(2));

  // The residual of the surface flow equation provides the mass flux from
  // subsurface to surface.
  Epetra_MultiVector& source = *S_next_->GetFieldData("surface_subsurface_flux",
          name_)->ViewComponent("cell",false);
  source = *g->SubVector(2)->Data()->ViewComponent("cell",false);

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->Functional(t_old, t_new, u_old->SubVector(0),
          u_new->SubVector(0), g->SubVector(0));

  // All surface to subsurface fluxes have been taken by the subsurface.
  g->SubVector(2)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  // Now that mass fluxes are done, do energy.
  // Evaluate the surface energy residual
  surf_energy_pk_->Functional(t_old, t_new, u_old->SubVector(3),
          u_new->SubVector(3), g->SubVector(3));

  // The residual of the surface energy equation provides the diffusive energy
  // flux from subsurface to surface.
  Epetra_MultiVector& esource =
      *S_next_->GetFieldData("surface_subsurface_energy_flux", name_)
      ->ViewComponent("cell",false);
  esource = *g->SubVector(3)->Data()->ViewComponent("cell",false);

  // Evaluate the subsurface energy residual.
  domain_energy_pk_->Functional(t_old, t_new, u_old->SubVector(1),
          u_new->SubVector(1), g->SubVector(1));

  // All energy fluxes have been taken by the subsurface.
  g->SubVector(3)->Data()->ViewComponent("cell",false)->PutScalar(0.);
}

// -- Apply preconditioner
int MPCPermafrost3::ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
        Teuchos::RCP<TreeVector> Pr) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon application:" << std::endl;

  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Residuals (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  r_ps");
    vnames.push_back("  r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(r->SubVector(2)->Data().ptr());
    vecs.push_back(r->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Residuals (subsurface):" << std::endl;
    vnames[0] = "  r_p";
    vnames[1] = "  r_T";
    vecs[0] = r->SubVector(0)->Data().ptr();
    vecs[1] = r->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  // make a new TreeVector that is just the subsurface values (by pointer).
  // -- note these const casts are necessary to create the new TreeVector, but
  //    since the TreeVector COULD be const (it is only used in a single method,
  //    in which it is const), const-correctness is not violated here.
  Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector());
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(0)));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(1)));

  Teuchos::RCP<TreeVector> domain_Pu_tv = Teuchos::rcp(new TreeVector());
  domain_Pu_tv->PushBack(Pr->SubVector(0));
  domain_Pu_tv->PushBack(Pr->SubVector(1));

  // call the operator's inverse
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying coupled subsurface operator." << std::endl;
  int ierr = lin_solver_->ApplyInverse(*domain_u_tv, *domain_Pu_tv);
  ierr = (ierr > 0) ? 0 : 1;

  if (S_next_->cycle() == 339) {
    Epetra_Vector vec(precon_->Schur()->Map());
    precon_->Schur()->ExtractDiagonalCopy(vec);
    vec.Print(std::cout);
  }

  
  // call EWC precon
  if (precon_type_ == PRECON_EWC) {
    // dump std correction to screen
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "PC_std * residuals (surface):" << std::endl;
      std::vector<std::string> vnames;
      vnames.push_back("  PC_std * r_ps");
      vnames.push_back("  PC_std * r_Ts");
      std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
      vecs.push_back(Pr->SubVector(2)->Data().ptr());
      vecs.push_back(Pr->SubVector(3)->Data().ptr());
      surf_db_->WriteVectors(vnames, vecs, true);

      *vo_->os() << "PC_std * residuals (subsurface):" << std::endl;
      vnames[0] = "  PC_std * r_p";
      vnames[1] = "  PC_std * r_T";
      vecs[0] = Pr->SubVector(0)->Data().ptr();
      vecs[1] = Pr->SubVector(1)->Data().ptr();
      domain_db_->WriteVectors(vnames, vecs, true);
    }

    // make sure we can back-calc face corrections that preserve residuals on faces
    Teuchos::RCP<TreeVector> res0 = Teuchos::rcp(new TreeVector(*domain_u_tv));
    res0->PutScalar(0.);
    Teuchos::RCP<TreeVector> Pu_std = Teuchos::rcp(new TreeVector(*domain_Pu_tv));
    *Pu_std = *domain_Pu_tv;

    // call EWC, which does Pu_p <-- Pu_p_std + dPu_p
    int ierr_ewc = sub_ewc_->ApplyPreconditioner(domain_u_tv, domain_Pu_tv);
    ierr += ierr_ewc;
    
    // calculate dPu_lambda from dPu_p
    Pu_std->Update(1.0, *domain_Pu_tv, -1.0);
    precon_->UpdateConsistentFaceCorrection(*res0, Pu_std.ptr());

    // update Pu_lambda <-- Pu_lambda_std + dPu_lambda
    domain_Pu_tv->SubVector(0)->Data()->ViewComponent("face",false)->Update(1.,
            *Pu_std->SubVector(0)->Data()->ViewComponent("face",false), 1.);
    domain_Pu_tv->SubVector(1)->Data()->ViewComponent("face",false)->Update(1.,
            *Pu_std->SubVector(1)->Data()->ViewComponent("face",false), 1.);
  }
  // end EWC PRECON

  
  // Copy subsurface face corrections to surface cell corrections
  CopySubsurfaceToSurface(*Pr->SubVector(0)->Data(),
                          Pr->SubVector(2)->Data().ptr());
  CopySubsurfaceToSurface(*Pr->SubVector(1)->Data(),
                          Pr->SubVector(3)->Data().ptr());

  // Derive surface face corrections.
  pc_surf_energy_->UpdateConsistentFaceCorrection(*r->SubVector(3)->Data(),
          Pr->SubVector(3)->Data().ptr());
  UpdateConsistentFaceCorrectionWater_(r.ptr(), Teuchos::null, Pr.ptr());

  // dump to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "PC * residuals (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  PC * r_ps");
    vnames.push_back("  PC * r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(Pr->SubVector(2)->Data().ptr());
    vecs.push_back(Pr->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "PC * residuals (subsurface):" << std::endl;
    vnames[0] = "  PC * r_p";
    vnames[1] = "  PC * r_T";
    vecs[0] = Pr->SubVector(0)->Data().ptr();
    vecs[1] = Pr->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }
  
  return ierr;
}

// -- Update the preconditioner.
void
MPCPermafrost3::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // Create the on-diagonal block PCs
  //  StrongMPC<PKPhysicalBDFBase>::UpdatePreconditioner(t,up,h);
  sub_pks_[2]->UpdatePreconditioner(t, up->SubVector(2), h);
  sub_pks_[3]->UpdatePreconditioner(t, up->SubVector(3), h);
  sub_pks_[0]->UpdatePreconditioner(t, up->SubVector(0), h);
  sub_pks_[1]->UpdatePreconditioner(t, up->SubVector(1), h);

  // Add the off-diagonal blocks.
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "temperature");
  Teuchos::RCP<const CompositeVector> dWCdT_domain =
      S_next_->GetFieldData("dwater_content_dtemperature");

  S_next_->GetFieldEvaluator("energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "pressure");
  Teuchos::RCP<const CompositeVector> dEdp_domain =
      S_next_->GetFieldData("denergy_dpressure");

  // write for debugging
  std::vector<std::string> vnames;
  vnames.push_back("  dwc_dT");
  vnames.push_back("  de_dp");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(dWCdT_domain.ptr());
  vecs.push_back(dEdp_domain.ptr());
  domain_db_->WriteVectors(vnames, vecs, false);

  // ALWAYS 0!
  // S_next_->GetFieldEvaluator("surface_water_content")
  //     ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_temperature");
  // Teuchos::RCP<const CompositeVector> dWCdT_surf =
  //     S_next_->GetFieldData("dsurface_water_content_dsurface_temperature");

  S_next_->GetFieldEvaluator("surface_energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_pressure");
  Teuchos::RCP<const CompositeVector> dEdp_surf =
      S_next_->GetFieldData("dsurface_energy_dsurface_pressure");

  // write for debugging
  vnames.resize(1);
  vecs.resize(1);
  vnames[0] = "  de_dp surf";
  vecs[0] = dEdp_surf.ptr();
  surf_db_->WriteVectors(vnames, vecs, false);

  // Scale by 1/h
  precon_->SetOffDiagonals(dWCdT_domain->ViewComponent("cell",false),
                           dEdp_domain->ViewComponent("cell",false),
                           Teuchos::null, // dWC_dT = 0
                           dEdp_surf->ViewComponent("cell",false),
                           1./h);

  // // update advective components
  // if (adv_flux_ == Teuchos::null) {
  //   Teuchos::RCP<Energy::EnergyBase> pk_as_energy = Teuchos::rcp_dynamic_cast<Energy::EnergyBase>(domain_energy_pk_);
  //   ASSERT(pk_as_energy != Teuchos::null);
  //   adv_flux_ = pk_as_energy->advection()->flux();
  // }

  // //     if (dynamic_mesh_) *pcAdv_ = *sub_pks_[0]->preconditioner();
  // Teuchos::RCP<const CompositeVector> rel_perm =
  //     S_next_->GetFieldData("numerical_rel_perm");
  // Teuchos::RCP<CompositeVector> rel_perm_times_enthalpy =
  //     Teuchos::rcp(new CompositeVector(*rel_perm));
  // *rel_perm_times_enthalpy = *rel_perm;
  // {
  //   Epetra_MultiVector& kr_f = *rel_perm_times_enthalpy->ViewComponent("face",false);
  //   const Epetra_MultiVector& enth_u = *adv_field_->ViewComponent("face",false);
  //   const Epetra_MultiVector& enth_c = *adv_field_->ViewComponent("cell",true);
  //   const Epetra_MultiVector& flux = *adv_flux_->ViewComponent("face",false);
  //   for (int f=0; f!=kr_f.MyLength(); ++f) {
  //     if (std::abs(flux[0][f]) > 1.e-12) {
  //       kr_f[0][f] *= enth_u[0][f] / std::abs(flux[0][f]);
  //     } else {
  //       AmanziMesh::Entity_ID_List cells;
  //       domain_mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  //       if (cells.size() == 1) {
  //         kr_f[0][f] *= enth_c[0][cells[0]];
  //       } else {
  //         kr_f[0][f] *= (enth_c[0][cells[0]] + enth_c[0][cells[1]])/2.;
  //       }
  //     }
  //   }
  // }
  // pcAdv_->CreateMFDstiffnessMatrices(rel_perm_times_enthalpy.ptr());
  
  // update ewc Precons if needed
  sub_ewc_->UpdatePreconditioner(t, up, h);
  surf_ewc_->UpdatePreconditioner(t, up, h);
}

// -- Modify the predictor.
bool
MPCPermafrost3::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  bool modified = false;

  // Make a new TreeVector that is just the subsurface values (by pointer).
  Teuchos::RCP<TreeVector> sub_u = Teuchos::rcp(new TreeVector());
  sub_u->PushBack(u->SubVector(0));
  sub_u->PushBack(u->SubVector(1));

  // Subsurface EWC, modifies cells
  modified |= sub_ewc_->ModifyPredictor(h,sub_u);

  // Calculate consistent faces
  modified |= domain_flow_pk_->ModifyPredictor(h, u0->SubVector(0), u->SubVector(0));
  modified |= domain_energy_pk_->ModifyPredictor(h, u0->SubVector(1), u->SubVector(1));

  // Copy consistent faces to surface
  if (modified) {
    S_next_->GetFieldEvaluator("surface_relative_permeability")->HasFieldChanged(S_next_.ptr(),name_);
    Teuchos::RCP<const CompositeVector> h_prev = S_inter_->GetFieldData("ponded_depth");
    MergeSubsurfaceAndSurfacePressure(*h_prev, u->SubVector(0)->Data().ptr(), u->SubVector(2)->Data().ptr());
    CopySubsurfaceToSurface(*u->SubVector(1)->Data(), u->SubVector(3)->Data().ptr());
  }

  // Hack surface faces
  bool newly_modified = false;
  newly_modified |= water_->ModifyPredictor_Heuristic(h, u);
  newly_modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
  newly_modified |= water_->ModifyPredictor_TempFromSource(h, u);
  modified |= newly_modified;

  // -- copy surf --> sub
  //  if (newly_modified) {
  CopySurfaceToSubsurface(*u->SubVector(2)->Data(), u->SubVector(0)->Data().ptr());
  CopySurfaceToSubsurface(*u->SubVector(3)->Data(), u->SubVector(1)->Data().ptr());
  //  }

  // Calculate consistent surface faces
  sub_pks_[2]->ChangedSolution();
  sub_pks_[3]->ChangedSolution();
  modified |= surf_flow_pk_->ModifyPredictor(h, u0->SubVector(2), u->SubVector(2));
  modified |= surf_energy_pk_->ModifyPredictor(h, u0->SubVector(3), u->SubVector(3));

  return modified;
}

// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCPermafrost3::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> r,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // dump NKAd correction to screen
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "NKA * PC * residuals (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  NKA*PC*r_ps");
    vnames.push_back("  NKA*PC*r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(du->SubVector(2)->Data().ptr());
    vecs.push_back(du->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "NKA * PC * residuals (subsurface):" << std::endl;
    vnames[0] = "  NKA*PC*r_p";
    vnames[1] = "  NKA*PC*r_T";
    vecs[0] = du->SubVector(0)->Data().ptr();
    vecs[1] = du->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  // modify correction using water approaches
  int n_modified = 0;
  n_modified += water_->ModifyCorrection_WaterFaceLimiter(h, r, u, du);
  double damping = water_->ModifyCorrection_WaterSpurtDamp(h, r, u, du);
  n_modified += water_->ModifyCorrection_WaterSpurtCap(h, r, u, du, damping);

  // -- accumulate globally
  int n_modified_l = n_modified;
  u->SubVector(0)->Data()->Comm().SumAll(&n_modified_l, &n_modified, 1);
  bool modified = (n_modified > 0) || (damping < 1.);

  // -- calculate consistent subsurface cells
  if (modified) {
    if (consistent_cells_) {
      // Derive subsurface cell corrections.
      pc_flow_->UpdateConsistentCellCorrection(
          *u->SubVector(0)->Data(),
          du->SubVector(0)->Data().ptr());
    }
  }

  // // modify correction on subsurface cells using EWC
  // if (precon_type_ == PRECON_EWC) {
  //     // make a new TreeVector that is just the subsurface values (by pointer).
  //     // -- note these const casts are necessary to create the new TreeVector, but
  //     //    since the TreeVector COULD be const (it is only used in a single method,
  //     //    in which it is const), const-correctness is not violated here.
  //     Teuchos::RCP<TreeVector> domain_res_tv = Teuchos::rcp(new TreeVector());
  //     domain_res_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(0)));
  //     domain_res_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(1)));
      
  //     Teuchos::RCP<TreeVector> domain_du_tv = Teuchos::rcp(new TreeVector());
  //     domain_du_tv->PushBack(du->SubVector(0));
  //     domain_du_tv->PushBack(du->SubVector(1));

  //     // make sure we can back-calc face corrections that preserve residuals on faces
  //     Teuchos::RCP<TreeVector> res0 = Teuchos::rcp(new TreeVector(*domain_res_tv));
  //     res0->PutScalar(0.);
  //     Teuchos::RCP<TreeVector> du_std = Teuchos::rcp(new TreeVector(*domain_du_tv));
  //     *du_std = *domain_du_tv;

  //     // call EWC, which does Pu_p <-- Pu_p_std + dPu_p
  //     sub_ewc_->ApplyPreconditioner(domain_res_tv, domain_du_tv);

  //     // calculate dPu_lambda from dPu_p
  //     du_std->Update(1.0, *domain_du_tv, -1.0);
  //     precon_->UpdateConsistentFaceCorrection(*res0, du_std.ptr());

  //     // update Pu_lambda <-- Pu_lambda_std + dPu_lambda
  //     domain_du_tv->SubVector(0)->Data()->ViewComponent("face",false)->Update(1.,
  //             *du_std->SubVector(0)->Data()->ViewComponent("face",false), 1.);
  //     domain_du_tv->SubVector(1)->Data()->ViewComponent("face",false)->Update(1.,
  //             *du_std->SubVector(1)->Data()->ViewComponent("face",false), 1.);
  //     modified = true;
  // }
  
  if (modified) {
    // Copy subsurface face corrections to surface cell corrections
    CopySubsurfaceToSurface(*du->SubVector(0)->Data(),
                            du->SubVector(2)->Data().ptr());

    // Derive surface face corrections.
    UpdateConsistentFaceCorrectionWater_(r.ptr(), u.ptr(), du.ptr());
  }

  // dump modified correction to screen
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Modified correction:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  Mod NKA*PC*r_ps");
    vnames.push_back("  Mod NKA*PC*r_Ts");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(du->SubVector(2)->Data().ptr());
    vecs.push_back(du->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Modified correction:" << std::endl;
    vnames[0] = "  Mod NKA*PC*r_p";
    vnames[1] = "  Mod NKA*PC*r_T";
    vecs[0] = du->SubVector(0)->Data().ptr();
    vecs[1] = du->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  return modified ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING :
      AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

void
MPCPermafrost3::UpdateConsistentFaceCorrectionWater_(const Teuchos::Ptr<const TreeVector>& r,
        const Teuchos::Ptr<const TreeVector>& u,
        const Teuchos::Ptr<TreeVector>& du) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // pull out corrections
  Teuchos::RCP<CompositeVector> surf_dp = du->SubVector(2)->Data();
  Epetra_MultiVector& surf_dp_c = *surf_dp->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> surf_dT = du->SubVector(3)->Data();
  Epetra_MultiVector& surf_dT_c = *surf_dT->ViewComponent("cell",false);

  // Calculate resulting delta h on the surface
  Teuchos::RCP<CompositeVector> surf_dh = Teuchos::rcp(new CompositeVector(*surf_dp));
  surf_dh->PutScalar(0.);

  // Grab the necessary components of the solution, in both CV and TV forms.
  Teuchos::RCP<CompositeVector> cv_p = S_next_->GetFieldData("surface_pressure", sub_pks_[2]->name());
  Teuchos::RCP<CompositeVector> cv_T = S_next_->GetFieldData("surface_temperature", sub_pks_[3]->name());

  Teuchos::RCP<const TreeVector> tv_p;
  Teuchos::RCP<const TreeVector> tv_T;
  if (u == Teuchos::null) {
    // This method is called by the ApplyPreconditioner(), which does not have access to u.
    Teuchos::RCP<TreeVector> tv_p_tmp = Teuchos::rcp(new TreeVector());
    tv_p_tmp->SetData(cv_p);
    tv_p = tv_p_tmp;

    Teuchos::RCP<TreeVector> tv_T_tmp = Teuchos::rcp(new TreeVector());
    tv_T_tmp->SetData(cv_T);
    tv_T = tv_T_tmp;
  } else {
    // This method is called by ModifyCorrection(), which does have access to u.
    tv_p = u->SubVector(2);
    tv_T = u->SubVector(3);
  }

  // Calculate/get old ponded depth, first ensuring PC didn't result in inadmissible solution
  if (sub_pks_[2]->IsAdmissible(tv_p) &&
      sub_pks_[3]->IsAdmissible(tv_T)) {
    S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
    *surf_dh->ViewComponent("cell",false) = *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false);

    // new ponded depth (again, checking admissibility)
    cv_p->ViewComponent("cell",false)->Update(-1., surf_dp_c, 1.);
    cv_T->ViewComponent("cell",false)->Update(-1., surf_dT_c, 1.);
    sub_pks_[2]->ChangedSolution();
    sub_pks_[3]->ChangedSolution();

    if (sub_pks_[2]->IsAdmissible(tv_p) &&
        sub_pks_[3]->IsAdmissible(tv_T)) {
      S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);

      // put delta ponded depth into surf_dh_cell
      surf_dh->ViewComponent("cell",false)
          ->Update(-1., *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false), 1.);

      // update delta faces
      pc_surf_flow_->UpdateConsistentFaceCorrection(*r->SubVector(2)->Data(), surf_dh.ptr());
      *surf_dp->ViewComponent("face",false) = *surf_dh->ViewComponent("face",false);

      // dump to screen
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Precon water correction." << std::endl;

        std::vector<std::string> vnames;
        vnames.push_back("pd_new");
        vnames.push_back("delta_pd");

        std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
        vecs.push_back(S_next_->GetFieldData("ponded_depth").ptr());
        vecs.push_back(surf_dh.ptr());
        surf_db_->WriteVectors(vnames, vecs, true);
      }
    }

    // revert solution so we don't break things
    S_next_->GetFieldData("surface_pressure",sub_pks_[2]->name())
        ->ViewComponent("cell",false)->Update(1., surf_dp_c, 1.);
    sub_pks_[2]->ChangedSolution();
    S_next_->GetFieldData("surface_temperature",sub_pks_[3]->name())
        ->ViewComponent("cell",false)->Update(1., surf_dT_c, 1.);
    sub_pks_[3]->ChangedSolution();
  }
}


int
MPCPermafrost3::ModifyCorrection_FrozenSurface_(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  const Epetra_MultiVector& res_surf_c = *res->SubVector(2)->Data()->ViewComponent("cell",false);
  const Epetra_MultiVector& u_surf_c = *u->SubVector(2)->Data()->ViewComponent("cell",false);

  Epetra_MultiVector& du_surf_c = *du->SubVector(2)->Data()->ViewComponent("cell",false);
  Epetra_MultiVector& du_domain_f = *du->SubVector(0)->Data()->ViewComponent("face",false);

  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
  S_next_->GetFieldEvaluator("surface_water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface_pressure");
  const Epetra_MultiVector& dWC_dp_surf_c =
      *S_next_->GetFieldData("dsurface_water_content_dsurface_pressure")
      ->ViewComponent("cell",false);

  int nmodified = 0;
  for (unsigned int sc=0; sc!=du_surf_c.MyLength(); ++sc) {
    if (std::abs(res_surf_c[0][sc]) > 0.) {
      nmodified++;

      du_surf_c[0][sc] = u_surf_c[0][sc] - (patm - res_surf_c[0][sc] / (dWC_dp_surf_c[0][sc] / h));
      AmanziMesh::Entity_ID f =
          surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
      du_domain_f[0][f] = du_surf_c[0][sc];
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  Source on frozen surface: res = " << res_surf_c[0][sc] << ", dWC_dp = " << dWC_dp_surf_c[0][sc]
                   << ", du = " << du_surf_c[0][sc] << std::endl;
    }
  }
  return nmodified;
}


} // namespace
