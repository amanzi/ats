
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"
#include "Operator_CellBndFace.hh"
#include "mpc_delegate_ewc_subsurface.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "permafrost_model.hh"
#include "surface_ice_model.hh"
#include "energy_base.hh"
#include "advection.hh"

#include "mpc_permafrost.hh"

namespace Amanzi {


MPCPermafrost::MPCPermafrost(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution) :
    PK(pk_tree, global_plist, S, solution),
    MPCSubsurface(pk_tree, global_plist, S, solution) {
  // tweak the sub-PK parameter lists
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");

  domain_subsurf_ = plist_->get<std::string>("domain name", "domain");
  if (domain_subsurf_ == "domain" || domain_subsurf_ == "") {
    domain_surf_ = plist_->get<std::string>("surface domain name", "surface");
  } else {
    domain_surf_ = plist_->get<std::string>("surface domain name", "surface_"+domain_subsurf_);
  }

  // propagate domain information down to delegates
  plist_->sublist("surface ewc delegate").set("domain name", domain_surf_);
  plist_->sublist("ewc delegate").set("domain name", domain_subsurf_);

  // exchange flux keys and evaluators
  mass_exchange_key_ = Keys::readKey(*plist_, domain_surf_, "mass exchange flux", "surface_subsurface_flux");
  energy_exchange_key_ = Keys::readKey(*plist_, domain_surf_, "energy exchange flux", "surface_subsurface_energy_flux");
  S->FEList().sublist(mass_exchange_key_).set("field evaluator type", "primary variable");
  S->FEList().sublist(energy_exchange_key_).set("field evaluator type", "primary variable");
}



void
MPCPermafrost::Setup(const Teuchos::Ptr<State>& S) {
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string> >("PKs order");
  // -- turn on coupling
  pks_list_->sublist(names[0]).set("coupled to surface via flux", true);
  pks_list_->sublist(names[1]).set("coupled to surface via flux", true);
  pks_list_->sublist(names[2]).set("coupled to subsurface via flux", true);
  pks_list_->sublist(names[3]).set("coupled to subsurface via flux", true);

  // -- ensure local ops are suface ops
  pks_list_->sublist(names[2]).sublist("diffusion preconditioner").set("surface operator", true);
  pks_list_->sublist(names[2]).sublist("accumulation preconditioner").set("surface operator", true);
  pks_list_->sublist(names[3]).sublist("diffusion preconditioner").set("surface operator", true);
  pks_list_->sublist(names[3]).sublist("advection preconditioner").set("surface operator", true);
  pks_list_->sublist(names[3]).sublist("accumulation preconditioner").set("surface operator", true);
  
  // grab the meshes
  surf_mesh_ = S->GetMesh(domain_surf_);
  domain_mesh_ = S->GetMesh(domain_subsurf_);
  
  // alias the PKs for easier reference
  domain_flow_pk_ = sub_pks_[0];
  domain_energy_pk_ = sub_pks_[1];
  surf_flow_pk_ = sub_pks_[2];
  surf_energy_pk_ = sub_pks_[3];

  
  // Create the dE_dp block, which will at least have a CELL-based diagonal
  // entry (from subsurface dE/dp) and a FACE-based diagonal entry (from
  // surface dE/dp), but the subsurface might only create a CELL-only matrix if
  // the other terms are supressed.  This can get removed/fixed once there is a
  // better way of creating/amalgamating ops into a single global operator.
  // For now this also means that we must have energy and flow using the same
  // discretization.
  Teuchos::ParameterList plist;
  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());

  std::string pk0_method = pks_list_->sublist(names[0]).sublist("diffusion").get<std::string>("discretization primary");
  std::string pk1_method = pks_list_->sublist(names[1]).sublist("diffusion").get<std::string>("discretization primary");
  if (pk0_method != pk1_method) {
    Errors::Message msg("MPC_Permafrost: for permafrost problems, due to issues in Jacobians, the flow and energy discretization methods must be the same.");
    Exceptions::amanzi_throw(msg);
  }
  
  if (pk0_method == "nlfv: bnd_faces" || pk0_method == "fv: bnd_faces") {
    cvs->SetMesh(domain_mesh_)->SetGhosted()
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1)
      ->AddComponent("cell", AmanziMesh::CELL, 1);  
  } else {
    cvs->SetMesh(domain_mesh_)->SetGhosted()
      ->AddComponent("face", AmanziMesh::FACE, 1)
      ->AddComponent("cell", AmanziMesh::CELL, 1);  
  }
  
  dE_dp_block_ = Teuchos::rcp(new Operators::Operator_FaceCell(cvs, plist));
  
  // call the subsurface setup, which calls the sub-pk's setups and sets up
  // the subsurface block operator
  MPCSubsurface::Setup(S);

  // require the coupling fields, claim ownership
  S->RequireField(mass_exchange_key_, name_)
      ->SetMesh(surf_mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::RCP<FieldEvaluator> fe = S->RequireFieldEvaluator(mass_exchange_key_);
  mass_exchange_pvfe_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe);
  AMANZI_ASSERT(mass_exchange_pvfe_.get());

  S->RequireField(energy_exchange_key_, name_)
      ->SetMesh(surf_mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  fe = S->RequireFieldEvaluator(energy_exchange_key_);
  energy_exchange_pvfe_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe);
  AMANZI_ASSERT(energy_exchange_pvfe_.get());
  
  if (precon_type_ != PRECON_NONE) {  
    // Add the (diagonal) surface blocks into the subsurface blocks.

    // For now we have just the basics, but this could get as complex as
    // MPCSubsurface with offdiagonal terms for surface advection, derivatives
    // of surface conductivity with respect to temperature, etc.

    // -- surface flow
    if (precon_type_ != PRECON_NO_FLOW_COUPLING) {
      Teuchos::RCP<Operators::Operator> surf_flow_pc = surf_flow_pk_->preconditioner();
      Teuchos::RCP<Operators::Operator> domain_flow_pc = domain_flow_pk_->preconditioner();
      for (Operators::Operator::op_iterator op = surf_flow_pc->begin();
           op != surf_flow_pc->end(); ++op) {
        domain_flow_pc->OpPushBack(*op);
      }
    }

    // -- surface energy
    Teuchos::RCP<Operators::Operator> surf_energy_pc = surf_energy_pk_->preconditioner();
    Teuchos::RCP<Operators::Operator> domain_energy_pc = domain_energy_pk_->preconditioner();
    for (Operators::Operator::op_iterator op = surf_energy_pc->begin();
         op != surf_energy_pc->end(); ++op) {
      domain_energy_pc->OpPushBack(*op);
    }

    if (precon_type_ != PRECON_BLOCK_DIAGONAL && precon_type_ != PRECON_NO_FLOW_COUPLING) {

      // Add off-diagonal blocks for the surface
      // -- derivatives of surface water content with respect to surface temperature
      // -- ALWAYS ZERO!
      // AMANZI_ASSERT(dWC_dT_block_ != Teuchos::null);
      // Teuchos::ParameterList dWC_dT_plist;
      // dWC_dT_plist.set("surface operator", true);
      // dWC_dT_plist.set("entity kind", "cell");
      // dWC_dT_surf_ = Teuchos::rcp(new Operators::PDE_Accumulation(dWC_dT_plist, dWC_dT_block_));

      // -- derivatives of surface energy with respect to surface pressure
      //    For the Operator, we have to create one with the surface mesh,
      //    then push the op into the full (subsurface) operator.
      AMANZI_ASSERT(dE_dp_block_ != Teuchos::null);
      Teuchos::ParameterList dE_dp_plist;
      dE_dp_plist.set("surface operator", true);
      dE_dp_plist.set("entity kind", "cell");
      dE_dp_surf_ = Teuchos::rcp(new Operators::PDE_Accumulation(dE_dp_plist, surf_mesh_));

      for (Operators::Operator::op_iterator op = dE_dp_surf_->global_operator()->begin();
           op != dE_dp_surf_->global_operator()->end(); ++op) {
        dE_dp_block_->OpPushBack(*op);
      }
    }

    // must now re-symbolic assemble the matrix to get the updated surface parts
    preconditioner_->set_inverse_parameters(plist_->sublist("preconditioner"));
  }
      
  // grab the debuggers
  domain_db_ = domain_flow_pk_->debugger();
  surf_db_ = surf_flow_pk_->debugger();

  // set up the Water delegate
  if (plist_->isSublist("water delegate")) {
    Teuchos::RCP<Teuchos::ParameterList> water_list = Teuchos::sublist(plist_, "water delegate");
    water_ = Teuchos::rcp(new MPCDelegateWater(water_list, domain_subsurf_));
    water_->set_indices(0,2,1,3);
    water_->set_db(surf_db_);
  }

}

void
MPCPermafrost::Initialize(const Teuchos::Ptr<State>& S) {
  // initialize coupling terms
  S->GetFieldData(mass_exchange_key_, name_)->PutScalar(0.);
  S->GetField(mass_exchange_key_, name_)->set_initialized();
  S->GetFieldData(energy_exchange_key_, name_)->PutScalar(0.);
  S->GetField(energy_exchange_key_, name_)->set_initialized();

  // Initialize all sub PKs.
  MPCSubsurface::Initialize(S);

  // ensure continuity of ICs... surface takes precedence if it was initialized
  if (S->GetField(Keys::getKey(domain_surf_, "pressure"))->initialized()) {
    CopySurfaceToSubsurface(*S->GetFieldData(Keys::getKey(domain_surf_,"pressure"), surf_flow_pk_->name()),
                            S->GetFieldData(Keys::getKey(domain_subsurf_,"pressure"), domain_flow_pk_->name()).ptr());
  } else {
    CopySubsurfaceToSurface(*S->GetFieldData(Keys::getKey(domain_subsurf_,"pressure"), domain_flow_pk_->name()),
                            S->GetFieldData(Keys::getKey(domain_surf_,"pressure"), surf_flow_pk_->name()).ptr());
  }
  if (S->GetField(Keys::getKey(domain_surf_, "temperature"))->initialized()) {
    CopySurfaceToSubsurface(*S->GetFieldData(Keys::getKey(domain_surf_,"temperature"), surf_energy_pk_->name()),
                            S->GetFieldData(Keys::getKey(domain_subsurf_,"temperature"), domain_energy_pk_->name()).ptr());
  } else {
    CopySubsurfaceToSurface(*S->GetFieldData(Keys::getKey(domain_subsurf_,"temperature"), domain_energy_pk_->name()),
                            S->GetFieldData(Keys::getKey(domain_surf_,"temperature"), surf_energy_pk_->name()).ptr());
  }
}


void
MPCPermafrost::set_states(const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<State>& S_inter,
                           const Teuchos::RCP<State>& S_next) {
  MPCSubsurface::set_states(S,S_inter,S_next);
  if (water_.get()) water_->set_states(S,S_inter,S_next);
}


// Compute the non-linear functional g = g(t,u,udot)
void
MPCPermafrost::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // propagate updated info into state
  Solution_to_State(*u_new, S_next_);

  // Evaluate the surface flow residual
  surf_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(2),
                            u_new->SubVector(2), g->SubVector(2));

  // The residual of the surface flow equation provides the mass flux from
  // subsurface to surface.
  Epetra_MultiVector& source = *S_next_->GetFieldData(mass_exchange_key_, name_)->ViewComponent("cell",false);
  source = *g->SubVector(2)->Data()->ViewComponent("cell",false);
  mass_exchange_pvfe_->SetFieldAsChanged(S_next_.ptr());

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  domain_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0),
          u_new->SubVector(0), g->SubVector(0));

  // All surface to subsurface fluxes have been taken by the subsurface.
  g->SubVector(2)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  // Now that mass fluxes are done, do energy.
  // Evaluate the surface energy residual
  surf_energy_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(3),
          u_new->SubVector(3), g->SubVector(3));

  // The residual of the surface energy equation provides the diffusive energy
  // flux from subsurface to surface.
  Epetra_MultiVector& esource =
      *S_next_->GetFieldData(energy_exchange_key_, name_)->ViewComponent("cell",false);
  esource = *g->SubVector(3)->Data()->ViewComponent("cell",false);
  energy_exchange_pvfe_->SetFieldAsChanged(S_next_.ptr());

  // Evaluate the subsurface energy residual.
  domain_energy_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(1),
          u_new->SubVector(1), g->SubVector(1));

  // All energy fluxes have been taken by the subsurface.
  g->SubVector(3)->Data()->ViewComponent("cell",false)->PutScalar(0.);
}

// -- Apply preconditioner
int MPCPermafrost::ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
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
  // since the TreeVector COULD be const (it is only used in a single method,
  // in which it is const), const-correctness is not violated here.  The
  // correct solution would be to have a TV constructor that took const
  // subvectors and made a const TV?
  Teuchos::RCP<TreeVector> domain_u_tv = Teuchos::rcp(new TreeVector());
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(0)));
  domain_u_tv->PushBack(Teuchos::rcp_const_cast<TreeVector>(r->SubVector(1)));

  Teuchos::RCP<TreeVector> domain_Pu_tv = Teuchos::rcp(new TreeVector());
  domain_Pu_tv->PushBack(Pr->SubVector(0));
  domain_Pu_tv->PushBack(Pr->SubVector(1));

  // call the operator's inverse
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon applying coupled subsurface operator." << std::endl;
  int ierr = preconditioner_->ApplyInverse(*domain_u_tv, *domain_Pu_tv);

  // rescale to Pa from MPa
  Pr->SubVector(0)->Data()->Scale(1.e6);

  // Copy subsurface face corrections to surface cell corrections
  CopySubsurfaceToSurface(*Pr->SubVector(0)->Data(),
                          Pr->SubVector(2)->Data().ptr());
  CopySubsurfaceToSurface(*Pr->SubVector(1)->Data(),
                          Pr->SubVector(3)->Data().ptr());

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
  }
  return (ierr > 0) ? 0 : 1;
}

// -- Update the preconditioner.
void
MPCPermafrost::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update the various components -- note it is important that subsurface are
  // done first (which is handled as they are listed first)
  MPCSubsurface::UpdatePreconditioner(t, up, h);
  
  // Add the surface off-diagonal blocks.
  // -- surface dWC/dT
  // -- ALWAYS ZERO
  // S_next_->GetFieldEvaluator("surface-water_content")
  //     ->HasFieldDerivativeChanged(S_next_.ptr(), name_, "surface-temperature");
  // Teuchos::RCP<const CompositeVector> dWCdT =
  //     S_next_->GetFieldData("dsurface-water_content_dsurface-temperature");
  // dWC_dT_surf_->AddAccumulationTerm(*dWCdT->ViewComponent("cell", false), h);
  

  if (precon_type_ != PRECON_NO_FLOW_COUPLING) {
    // -- surface dE_dp
    S_next_->GetFieldEvaluator(Keys::getKey(domain_surf_,"energy"))
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, Keys::getKey(domain_surf_,"pressure"));
    Teuchos::RCP<const CompositeVector> dEdp =
        S_next_->GetFieldData(Keys::getDerivKey(Keys::getKey(domain_surf_,"energy"), Keys::getKey(domain_surf_,"pressure")));
    dE_dp_surf_->AddAccumulationTerm(*dEdp, h, "cell", false);
  
    // write for debugging
    std::vector<std::string> vnames;
    //  vnames.push_back("  dwc_dT");
    vnames.push_back("  de_dp");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    //  vecs.push_back(dWCdT.ptr());
    vecs.push_back(dEdp.ptr());
    surf_db_->WriteVectors(vnames, vecs, false);
  }

  // assemble
  // -- scale the pressure dofs
  double scaling = 1.e6; // dWC/dp_Pa * (Pa / MPa) --> dWC/dp_MPa
  sub_pks_[0]->preconditioner()->Rescale(scaling);
  dE_dp_block_->Rescale(scaling);
  

  if (dump_) {
    std::stringstream filename;
    filename << "FullyCoupled_PC_" << S_next_->cycle() << "_" << update_pcs_ << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename.str().c_str(), *preconditioner_->A());
  }
  
  // update ewc Precons if needed
  //  surf_ewc_->UpdatePreconditioner(t, up, h);
}

// -- Modify the predictor.
bool
MPCPermafrost::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  bool modified = false;

  // HACK to allow for predictor use in subcycling, but then trash the history
  // if operator splitting coupler has overwritten our OLD time's value
  S_next_->GetFieldEvaluator(Keys::getKey(domain_subsurf_, "water_content"))
    ->HasFieldChanged(S_next_.ptr(), name_);
  if (S_inter_->GetFieldEvaluator(Keys::getKey(domain_subsurf_, "water_content"))
      ->HasFieldChanged(S_inter_.ptr(), name_)) {
    *u = *u0;
    ChangedSolution();
    S_next_->GetFieldEvaluator(Keys::getKey(domain_subsurf_, "water_content"))
      ->HasFieldChanged(S_next_.ptr(), name_);
    return false; // intentionally lieing -- true here triggers another call of ChangedSolution() which we want to avoid
  }
  
  // write predictor
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Extrapolated Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_extrap");
    vnames.push_back("  Ts_extrap");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Extrapolated Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_extrap";
    vnames[1] = "  T_extrap";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }

  // Make a new TreeVector that is just the subsurface values (by pointer).
  Teuchos::RCP<TreeVector> sub_u = Teuchos::rcp(new TreeVector());
  sub_u->PushBack(u->SubVector(0));
  sub_u->PushBack(u->SubVector(1));

  // Subsurface EWC, modifies cells
  modified |= ewc_->ModifyPredictor(h,sub_u);

  // write predictor
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "EWC Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_ewc");
    vnames.push_back("  Ts_ewc");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "EWC Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_ewc";
    vnames[1] = "  T_ewc";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }
 
  // Calculate consistent faces
  modified |= domain_flow_pk_->ModifyPredictor(h, u0->SubVector(0), u->SubVector(0));
  modified |= domain_energy_pk_->ModifyPredictor(h, u0->SubVector(1), u->SubVector(1));

  // write predictor
  if (modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "EWC/Consistent Face Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_extrap");
    vnames.push_back("  Ts_extrap");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "EWC/Consistent Face Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_extrap";
    vnames[1] = "  T_extrap";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }
  
  // Copy consistent faces to surface
  if (modified) {
    //S_next_->GetFieldEvaluator(Keys::getKey(domain_surf_,"relative_permeability"))->HasFieldChanged(S_next_.ptr(),name_);
    Teuchos::RCP<const CompositeVector> h_prev = S_inter_->GetFieldData(Keys::getKey(domain_surf_,"ponded_depth"));

    MergeSubsurfaceAndSurfacePressure(*h_prev, u->SubVector(0)->Data().ptr(), u->SubVector(2)->Data().ptr());
    CopySubsurfaceToSurface(*u->SubVector(1)->Data(), u->SubVector(3)->Data().ptr());

  }

  // Hack surface faces
  bool newly_modified = false;
  if (water_.get()) {
    newly_modified |= water_->ModifyPredictor_Heuristic(h, u);
    newly_modified |= water_->ModifyPredictor_WaterSpurtDamp(h, u);
    newly_modified |= water_->ModifyPredictor_TempFromSource(h, u);
    modified |= newly_modified;
  }

  // write predictor
  if (newly_modified && vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Spurt Fixed Prediction (surface):" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  ps_spurt");
    vnames.push_back("  Ts_spurt");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(2)->Data().ptr());
    vecs.push_back(u->SubVector(3)->Data().ptr());
    surf_db_->WriteVectors(vnames, vecs, true);

    *vo_->os() << "Spurt Fixed Prediction (subsurface):" << std::endl;
    vnames[0] = "  p_spurt";
    vnames[1] = "  T_spurt";
    vecs[0] = u->SubVector(0)->Data().ptr();
    vecs[1] = u->SubVector(1)->Data().ptr();
    domain_db_->WriteVectors(vnames, vecs, true);
  }
  

  // -- copy surf --> sub
  //  if (newly_modified) {
  CopySurfaceToSubsurface(*u->SubVector(2)->Data(), u->SubVector(0)->Data().ptr());
  CopySurfaceToSubsurface(*u->SubVector(3)->Data(), u->SubVector(1)->Data().ptr());
  //  }

  // Calculate consistent surface faces
  surf_flow_pk_->ChangedSolution();
  surf_energy_pk_->ChangedSolution();
  modified |= surf_flow_pk_->ModifyPredictor(h, u0->SubVector(2), u->SubVector(2));
  modified |= surf_energy_pk_->ModifyPredictor(h, u0->SubVector(3), u->SubVector(3));
  
  return modified;
}

// -- Modify the correction.
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
MPCPermafrost::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> r,
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

  // apply PK modifications
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult pk_modified =   
      StrongMPC<PK_PhysicalBDF_Default>::ModifyCorrection(h,r,u,du);
  if (pk_modified) {
    CopySurfaceToSubsurface(*du->SubVector(2)->Data(),
                            du->SubVector(0)->Data().ptr());
    CopySurfaceToSubsurface(*du->SubVector(3)->Data(),
                            du->SubVector(1)->Data().ptr());
  }
    
  // modify correction using water approaches
  int n_modified = 0;
  double damping = 1;
  if (water_.get()) {
    damping = water_->ModifyCorrection_SaturatedSpurtDamp(h, r, u, du);
    n_modified += water_->ModifyCorrection_SaturatedSpurtCap(h, r, u, du, damping);

    double damping_surf = water_->ModifyCorrection_WaterSpurtDamp(h, r, u, du);
    n_modified += water_->ModifyCorrection_WaterSpurtCap(h, r, u, du, damping_surf);

    // -- total damping
    damping = damping * damping_surf;
  
    // -- accumulate globally
    int n_modified_l = n_modified;
    u->SubVector(0)->Data()->Comm()->SumAll(&n_modified_l, &n_modified, 1);
  }
  bool modified = (n_modified > 0) || (damping < 1.);

  if (modified) {
    // Copy subsurface face corrections to surface cell corrections
    CopySubsurfaceToSurface(*du->SubVector(0)->Data(),
                            du->SubVector(2)->Data().ptr());
  }

  // dump modified correction to screen
  if ((modified || pk_modified) && vo_->os_OK(Teuchos::VERB_HIGH)) {
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

  if (modified) {
    // disallow backtracking which takes us back under patm
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING;
  } else {
    return pk_modified;
  }
}


} // namespace
