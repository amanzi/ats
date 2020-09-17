/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for coupling energy and water in the subsurface,
with freezing.

------------------------------------------------------------------------- */
#include "EpetraExt_RowMatrixOut.h"

#include "MultiplicativeEvaluator.hh"
#include "TreeOperator.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Advection.hh"
#include "PDE_Accumulation.hh"
#include "Operator.hh"
#include "upwind_total_flux.hh"
#include "upwind_arithmetic_mean.hh"

#include "permafrost_model.hh"
#include "liquid_ice_model.hh"
#include "richards.hh"

#include "mpc_delegate_ewc_subsurface.hh"
#include "mpc_surface.hh"

#define DEBUG_FLAG 1

namespace Amanzi {

// -- Initialize owned (dependent) variables.
void MPCSurface::Setup(const Teuchos::Ptr<State>& S) {
  // set up keys

  ddivq_dT_key_ = "ddivq_dT";
  ddivKgT_dp_key_ = "ddivKgT_dp";

  auto pk_order = plist_->get< Teuchos::Array<std::string> >("PKs order");
  domain_name_ = plist_->get<std::string>("domain name");
 
  temp_key_ = Keys::readKey(*plist_, domain_name_, "temperature", "temperature");
  pres_key_ = Keys::readKey(*plist_, domain_name_, "pressure", "pressure");
  e_key_ = Keys::readKey(*plist_, domain_name_, "energy", "energy");
  wc_key_ = Keys::readKey(*plist_, domain_name_, "water content", "water_content");

  kr_key_ = Keys::readKey(*plist_, domain_name_, "overland conductivity", "overland_conductivity");
  
  // set up the sub-pks
  StrongMPC<PK_PhysicalBDF_Default>::Setup(S);
  mesh_ = S->GetMesh(domain_name_);

  // set up debugger
  db_ = sub_pks_[0]->debugger();

  // Get the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Operator> pcA = sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Operator> pcB = sub_pks_[1]->preconditioner();

  // Create the combined operator
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcA->DomainMap()))));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcB->DomainMap()))));
 
  preconditioner_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  preconditioner_->SetOperatorBlock(0, 0, pcA);
  preconditioner_->SetOperatorBlock(1, 1, pcB);
  
  // select the method used for preconditioning
  std::string precon_string = plist_->get<std::string>("preconditioner type",
                                                       "picard");
  if (precon_string == "none") {
    precon_type_ = PRECON_NONE;
  } else if (precon_string == "block diagonal") {
    precon_type_ = PRECON_BLOCK_DIAGONAL;
  } else if (precon_string == "no flow coupling") {
    precon_type_ = PRECON_NO_FLOW_COUPLING;
  } else if (precon_string == "picard") {
    precon_type_ = PRECON_PICARD;
  } else if (precon_string == "ewc") {
    AMANZI_ASSERT(0);
    precon_type_ = PRECON_EWC;
  } else if (precon_string == "smart ewc") {
    AMANZI_ASSERT(0);
    precon_type_ = PRECON_EWC;
  } else {
    Errors::Message message(std::string("Invalid preconditioner type ")+precon_string);
    Exceptions::amanzi_throw(message);
  }

  // create offdiagonal blocks
  if (precon_type_ != PRECON_NONE && precon_type_ != PRECON_BLOCK_DIAGONAL) {
    // Create the block for derivatives of mass conservation with respect to temperature
    // -- derivatives of kr with respect to temperature
    if (precon_type_ != PRECON_NO_FLOW_COUPLING && 
        !plist_->get<bool>("supress Jacobian terms: d div q / dT", false)) {
      // set up the operator
      Teuchos::ParameterList divq_plist(pks_list_->sublist(pk_order[0]).sublist("diffusion preconditioner"));
      divq_plist.set("Newton correction", "true Jacobian");
      divq_plist.set("exclude primary terms", true);
      Operators::PDE_DiffusionFactory opfactory;
      ddivq_dT_ = opfactory.Create(divq_plist, mesh_);
      dWC_dT_block_ = ddivq_dT_->global_operator();
    }

    // -- derivatives of water content with respect to temperature
    if (dWC_dT_block_ == Teuchos::null) {
      dWC_dT_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
      dWC_dT_block_ = dWC_dT_->global_operator();
    } else {
      dWC_dT_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, dWC_dT_block_));
    }

    // // -- derivatives of advection term
    // if (precon_type_ != PRECON_NO_FLOW_COUPLING && 
    //     !plist_->get<bool>("supress Jacobian terms: div hq / dp,T", false)) {
    //   // derivative with respect to pressure
    //   Teuchos::ParameterList divhq_dp_plist(pks_list_->sublist(pk_order[0]).sublist("diffusion preconditioner"));

    //   divhq_dp_plist.set("Newton correction", "true Jacobian");

    //   Operators::PDE_DiffusionFactory opfactory;
    //   if (dE_dp_block_ == Teuchos::null) {
    //     ddivhq_dp_ = opfactory.Create(divhq_dp_plist, mesh_);
    //     dE_dp_block_ = ddivhq_dp_->global_operator();
    //   } else {
    //     ddivhq_dp_ = opfactory.Create(divhq_dp_plist, dE_dp_block_);
    //   }

    //   // derivative with respect to temperature
    //   Teuchos::ParameterList divhq_dT_plist(pks_list_->sublist(pk_order[0]).sublist("diffusion preconditioner"));
    //   divhq_dT_plist.set("exclude primary terms", true);

    //   divhq_dT_plist.set("Newton correction", "true Jacobian");
    //   ddivhq_dT_ = opfactory.Create(divhq_dT_plist, pcB);

    //   // need a field, evaluator, and upwinding for h * kr * rho/mu
    //   // -- first the evaluator
    //   Teuchos::ParameterList hkr_eval_list;
    //   hkr_eval_list.set("evaluator name", hkr_key_);
    //   Teuchos::Array<std::string> deps(2);
    //   deps[0] = enth_key_; deps[1] = kr_key_;
    //   hkr_eval_list.set("evaluator dependencies", deps);
    //   Teuchos::RCP<FieldEvaluator> hkr_eval =
    //       Teuchos::rcp(new Relations::MultiplicativeEvaluator(hkr_eval_list));

    //   // -- now the field
    //   names2[1] = "boundary_face";
    //   locations2[1] = AmanziMesh::BOUNDARY_FACE;
    //   S->RequireField(hkr_key_)->SetMesh(mesh_)->SetGhosted()
    //       ->AddComponents(names2, locations2, num_dofs2);
    //   S->SetFieldEvaluator(hkr_key_, hkr_eval);

    //   // locations2[1] = AmanziMesh::FACE;
    //   // names2[1] = "face";
    //   S->RequireField(uw_hkr_key_, name_)
    //     ->SetMesh(mesh_)->SetGhosted()
    //     ->SetComponent("face", AmanziMesh::FACE, 1);
    //       // ->SetComponents(names2, locations2, num_dofs2);
    //   S->GetField(uw_hkr_key_,name_)->set_io_vis(false);

    //   std::string method_name = pks_list_->sublist(pk_order[0])
    //       .get<std::string>("relative permeability method", "upwind with gravity");
    //   if (method_name != "upwind with Darcy flux") {
    //     Errors::Message msg;
    //     msg << "Subsurface coupler with advective Jacobian terms only supports a Richards upwind scheme of "
    //         << "\"upwind with Darcy flux\", but the method \"" << method_name << "\" was requested.";
    //     Exceptions::amanzi_throw(msg);
    //   }
    //   upwinding_hkr_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
    //           hkr_key_, uw_hkr_key_, mass_flux_dir_key_, 1.e-8));
      
    //   if (!is_fv_) {
    //     // -- and the upwinded field

    //     locations2[1] = AmanziMesh::FACE;
    //     names2[1] = "face";

    //     S->RequireField(Keys::getDerivKey(uw_hkr_key_, pres_key_), name_)
    //         ->SetMesh(mesh_)->SetGhosted()
    //         ->SetComponent("face", AmanziMesh::FACE, 1);
    //     S->GetField(Keys::getDerivKey(uw_hkr_key_, pres_key_),name_)
    //         ->set_io_vis(false);
    //     S->RequireField(Keys::getDerivKey(uw_hkr_key_, temp_key_), name_)
    //         ->SetMesh(mesh_)->SetGhosted()
    //         ->SetComponent("face", AmanziMesh::FACE, 1);
    //     S->GetField(Keys::getDerivKey(uw_hkr_key_, temp_key_), name_)
    //         ->set_io_vis(false);


    //     // -- and the upwinding
    //     upwinding_dhkr_dp_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
    //             Keys::getDerivKey(hkr_key_, pres_key_),
    //             Keys::getDerivKey(uw_hkr_key_, pres_key_),
    //             mass_flux_dir_key_, 1.e-8));
    //     upwinding_dhkr_dT_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
    //             Keys::getDerivKey(hkr_key_, temp_key_),
    //             Keys::getDerivKey(uw_hkr_key_, temp_key_),
    //             mass_flux_dir_key_, 1.e-8));
    //   }
    // }

    // // S->RequireField(ddivq_dT_key_, name_)->SetMesh(mesh_)->SetGhosted()
    // //   ->SetComponent("cell",AmanziMesh::CELL, 1);
    // // S->GetField(ddivq_dT_key_, name_)->set_io_vis(true);
    // // S->RequireField(ddivKgT_dp_key_, name_)->SetMesh(mesh_)->SetGhosted()
    // //   ->SetComponent("cell",AmanziMesh::CELL, 1);
    // // S->GetField(ddivKgT_dp_key_, name_)->set_io_vis(true);


    // -- derivatives of energy with respect to pressure
    if (dE_dp_block_ == Teuchos::null) {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
      dE_dp_block_ = dE_dp_->global_operator();
    } else {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, dE_dp_block_));
    }

    AMANZI_ASSERT(dWC_dT_block_ != Teuchos::null);
    AMANZI_ASSERT(dE_dp_block_ != Teuchos::null);
    preconditioner_->SetOperatorBlock(0, 1, dWC_dT_block_);
    preconditioner_->SetOperatorBlock(1, 0, dE_dp_block_);
  }

  // set up sparsity structure
  preconditioner_->set_inverse_parameters(plist_->sublist("preconditioner"));
}

void MPCSurface::Initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC<PK_PhysicalBDF_Default>::Initialize(S);

  // if (precon_type_ != PRECON_NONE && precon_type_ != PRECON_BLOCK_DIAGONAL) {
  //   if (S->HasField(ddivq_dT_key_)){
  //     S->GetFieldData(ddivq_dT_key_, name_)->PutScalar(0.0);
  //     S->GetField(ddivq_dT_key_, name_)->set_initialized();
  //   }
  // }

  if (ddivq_dT_ != Teuchos::null) {
    ddivq_dT_->SetBCs(sub_pks_[0]->BCs(), sub_pks_[1]->BCs());
    ddivq_dT_->SetTensorCoefficient(Teuchos::null);
  }

  // if (ddivhq_dp_ != Teuchos::null) {

  //   S->GetFieldData(uw_hkr_key_, name_)->PutScalar(1.);
  //   S->GetField(uw_hkr_key_, name_)->set_initialized();

  //   if (!is_fv_) {
  //     S->GetFieldData(Keys::getDerivKey(uw_hkr_key_, pres_key_), name_)->PutScalar(0.);
  //     S->GetField(Keys::getDerivKey(uw_hkr_key_, pres_key_), name_)->set_initialized();
  //     S->GetFieldData(Keys::getDerivKey(uw_hkr_key_, temp_key_), name_)->PutScalar(0.);
  //     S->GetField(Keys::getDerivKey(uw_hkr_key_, temp_key_), name_)->set_initialized();
  //   }

  //   Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
  //   AmanziGeometry::Point g(3);
  //   g[0] = (*gvec)[0]; g[1] = (*gvec)[1]; g[2] = (*gvec)[2];
  //   ddivhq_dp_->SetGravity(g);    
  //   ddivhq_dp_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[0]->BCs());
  //   ASSERT(richards_pk_ != Teuchos::null);
  //   ddivhq_dp_->SetTensorCoefficient(richards_pk_->K_);

  //   ddivhq_dT_->SetGravity(g);    
  //   ddivhq_dT_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[1]->BCs());
  //   ASSERT(richards_pk_ != Teuchos::null);
  //   ddivhq_dT_->SetTensorCoefficient(richards_pk_->K_);
  // }  

}


// updates the preconditioner
void MPCSurface::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (precon_type_ == PRECON_NONE) {
    // nothing to do
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    StrongMPC::UpdatePreconditioner(t,up,h);
  } else if (precon_type_ == PRECON_PICARD || precon_type_ == PRECON_EWC) {
    StrongMPC::UpdatePreconditioner(t,up,h);

    // Update operators for off-diagonals
    dWC_dT_block_->Init();
    dE_dp_block_->Init();
    
    // dWC / dT block
    // -- dkr/dT
    if (ddivq_dT_ != Teuchos::null) {
      // -- update and upwind d kr / dT
      S_next_->GetFieldEvaluator(kr_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, temp_key_);
      Teuchos::RCP<const CompositeVector> dkrdT;
      dkrdT = S_next_->GetFieldData(Keys::getDerivKey(kr_key_, temp_key_));

      // form the operator
      ddivq_dT_->SetScalarCoefficient(Teuchos::null, dkrdT);
      ddivq_dT_->UpdateMatrices(Teuchos::null,
              up->SubVector(0)->Data().ptr());              
      ddivq_dT_->UpdateMatricesNewtonCorrection(Teuchos::null,
              up->SubVector(0)->Data().ptr());              
      ddivq_dT_->ApplyBCs(false, true, false);
    }

    // -- dWC/dT diagonal term
    if (dWC_dT_ != Teuchos::null &&
        S_next_->GetFieldEvaluator(wc_key_)->IsDependency(S_next_.ptr(), temp_key_)) {
      S_next_->GetFieldEvaluator(wc_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, temp_key_);
      Teuchos::RCP<const CompositeVector> dWC_dT =
          S_next_->GetFieldData(Keys::getDerivKey(wc_key_, temp_key_));
      dWC_dT_->AddAccumulationTerm(*dWC_dT, h, "cell", false);

      // write for debugging
      db_->WriteVector("  dWC_dT", dWC_dT.ptr(), false);

    }

    // // -- d adv / dp   This one is a bit more complicated...
    // if (ddivhq_dp_ != Teuchos::null) {
    //   // Update and upwind enthalpy * kr * rho/mu
    //   // -- update values
    //   S_next_->GetFieldEvaluator(hkr_key_)
    //       ->HasFieldChanged(S_next_.ptr(), name_);
    //   S_next_->GetFieldEvaluator(hkr_key_)
    //       ->HasFieldDerivativeChanged(S_next_.ptr(), name_, pres_key_);
    //   S_next_->GetFieldEvaluator(hkr_key_)
    //       ->HasFieldDerivativeChanged(S_next_.ptr(), name_, temp_key_);

    //   Teuchos::RCP<const CompositeVector> denth_kr_dp_uw;
    //   Teuchos::RCP<const CompositeVector> denth_kr_dT_uw;

    //   Teuchos::RCP<const CompositeVector> enth_kr =
    //     S_next_->GetFieldData(hkr_key_);
    //   Teuchos::RCP<CompositeVector> enth_kr_uw =
    //     S_next_->GetFieldData(uw_hkr_key_, name_);
    //   enth_kr_uw->PutScalar(0.);

    //   enth_kr_uw->ViewComponent("face",false)
    //       ->Export(*enth_kr->ViewComponent("boundary_face",false),
    //                mesh_->exterior_face_importer(), Insert);
    //   upwinding_hkr_->Update(S_next_.ptr(), db_.ptr());

    //   ASSERT(richards_pk_ != Teuchos::null);
    //   if (richards_pk_->clobber_surf_kr_) {
    //     // -- stick zeros in the boundary faces
    //     Epetra_MultiVector enth_kr_bf(*enth_kr->ViewComponent("boundary_face",false));
    //     enth_kr_bf.PutScalar(0.0);
    //     enth_kr_uw->ViewComponent("face",false)->Export(enth_kr_bf,
    //             mesh_->exterior_face_importer(), Insert);
    //   }
      
    //   if (is_fv_) {
    //     denth_kr_dp_uw = S_next_->GetFieldData(Keys::getDerivKey(hkr_key_, pres_key_));
    //     denth_kr_dT_uw = S_next_->GetFieldData(Keys::getDerivKey(hkr_key_, temp_key_));
    //   } else {

    //     Teuchos::RCP<const CompositeVector> denth_kr_dp =
    //         S_next_->GetFieldData(Keys::getDerivKey(hkr_key_, pres_key_));
    //     Teuchos::RCP<const CompositeVector> denth_kr_dT =
    //         S_next_->GetFieldData(Keys::getDerivKey(hkr_key_, temp_key_));

    //     // -- zero target data (may be unnecessary?)
    //     Teuchos::RCP<CompositeVector> denth_kr_dp_uw_nc =
    //         S_next_->GetFieldData(Keys::getDerivKey(uw_hkr_key_, pres_key_), name_);
    //     denth_kr_dp_uw_nc->PutScalar(0.);
    //     Teuchos::RCP<CompositeVector> denth_kr_dT_uw_nc =
    //         S_next_->GetFieldData(Keys::getDerivKey(uw_hkr_key_, temp_key_), name_);
    //     denth_kr_dT_uw_nc->PutScalar(0.);

    //     // -- copy boundary faces into upwinded vector
    //     denth_kr_dp_uw_nc->ViewComponent("face",false)
    //         ->Export(*denth_kr_dp->ViewComponent("boundary_face",false),
    //                  mesh_->exterior_face_importer(), Insert);
    //     denth_kr_dT_uw_nc->ViewComponent("face",false)
    //         ->Export(*denth_kr_dT->ViewComponent("boundary_face",false),
    //                  mesh_->exterior_face_importer(), Insert);

    //     // -- upwind      
    //     upwinding_dhkr_dp_->Update(S_next_.ptr(), db_.ptr());
    //     upwinding_dhkr_dT_->Update(S_next_.ptr(), db_.ptr());

    //     // -- clobber
    //     ASSERT(richards_pk_ != Teuchos::null);
    //     if (richards_pk_->clobber_surf_kr_) {
    //       // -- stick zeros in the boundary faces
    //       Epetra_MultiVector enth_kr_bf(*enth_kr->ViewComponent("boundary_face",false));
    //       enth_kr_bf.PutScalar(0.0);
    //       denth_kr_dp_uw_nc->ViewComponent("face",false)->Export(enth_kr_bf,
    //               mesh_->exterior_face_importer(), Insert);
    //       denth_kr_dT_uw_nc->ViewComponent("face",false)->Export(enth_kr_bf,
    //               mesh_->exterior_face_importer(), Insert);
    //     }

    //     denth_kr_dp_uw =
    //         S_next_->GetFieldData(Keys::getDerivKey(uw_hkr_key_, pres_key_));
    //     denth_kr_dT_uw =
    //         S_next_->GetFieldData(Keys::getDerivKey(uw_hkr_key_, temp_key_));

    //   }
      
    //   Teuchos::RCP<const CompositeVector> flux = S_next_->GetFieldData(mass_flux_key_);
    //   Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(rho_key_);

    //   // form the operator: pressure component
    //   ddivhq_dp_->SetDensity(rho);
    //   ddivhq_dp_->SetScalarCoefficient(enth_kr_uw, denth_kr_dp_uw);
    //   // -- update the local matrices, div h * kr grad
    //   ddivhq_dp_->UpdateMatrices(Teuchos::null, Teuchos::null);
    //   // -- determine the advective fluxes, q_a = h * kr grad p
    //   CompositeVector adv_flux(*flux, INIT_MODE_ZERO);
    //   Teuchos::Ptr<CompositeVector> adv_flux_ptr(&adv_flux);
    //   ddivhq_dp_->UpdateFlux(up->SubVector(0)->Data().ptr(), adv_flux_ptr);
    //   // -- add in components div (d h*kr / dp) grad q_a / (h*kr)
    //   ddivhq_dp_->UpdateMatricesNewtonCorrection(adv_flux_ptr, up->SubVector(0)->Data().ptr());
    //   ddivhq_dp_->ApplyBCs(false, true, false);

    //   // form the operator: temperature component
    //   ddivhq_dT_->SetDensity(rho);
    //   ddivhq_dT_->SetScalarCoefficient(enth_kr_uw, denth_kr_dT_uw);
    //   // -- add in components div (d h*kr / dp) grad q_a / (h*kr)
    //   ddivhq_dT_->UpdateMatrices(adv_flux_ptr, up->SubVector(0)->Data().ptr());
    //   ddivhq_dT_->UpdateMatricesNewtonCorrection(adv_flux_ptr, up->SubVector(0)->Data().ptr());
    //   ddivhq_dT_->ApplyBCs(false, true, false);

    // }

    // -- dE/dp diagonal term
    S_next_->GetFieldEvaluator(e_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, pres_key_);
    Teuchos::RCP<const CompositeVector> dE_dp =
      S_next_->GetFieldData(Keys::getDerivKey(e_key_, pres_key_));
    dE_dp_->AddAccumulationTerm(*dE_dp, h, "cell", false);

    // write for debugging
    db_->WriteVector("  de_dp", dE_dp.ptr(), false);

    // finally assemble the full system, dump if requested, and form the inverse
  }
  update_pcs_++;
}


// -----------------------------------------------------------------------------
// Wrapper to call the requested preconditioner.
// -----------------------------------------------------------------------------
int MPCSurface::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon application:" << std::endl;

  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Residuals:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  r_p"); vnames.push_back("  r_T"); 
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(0)->Data().ptr()); 
    vecs.push_back(u->SubVector(1)->Data().ptr()); 
    db_->WriteVectors(vnames, vecs, true);
  }
  
  int ierr = 0;
  if (precon_type_ == PRECON_NONE) {
    *Pu = *u;
    ierr = 1;
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    ierr = StrongMPC::ApplyPreconditioner(u,Pu);
  } else if (precon_type_ == PRECON_PICARD) {
    ierr = preconditioner_->ApplyInverse(*u, *Pu);
  } else if (precon_type_ == PRECON_EWC) {
    ierr = preconditioner_->ApplyInverse(*u, *Pu);
  }

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "PC * residuals:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  PC*r_p"); vnames.push_back("  PC*r_T"); 
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(Pu->SubVector(0)->Data().ptr()); 
    vecs.push_back(Pu->SubVector(1)->Data().ptr()); 
    db_->WriteVectors(vnames, vecs, true);
  }
  
  return (ierr > 0) ? 0 : 1;
}



} // namespace
