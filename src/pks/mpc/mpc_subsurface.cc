/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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
#include "OperatorDiffusionFactory.hh"
#include "OperatorAdvection.hh"
#include "OperatorAccumulation.hh"
#include "Operator.hh"
#include "upwind_total_flux.hh"
#include "upwind_arithmetic_mean.hh"

#include "permafrost_model.hh"
#include "richards.hh"

#include "mpc_delegate_ewc_subsurface.hh"
#include "mpc_subsurface.hh"

#define DEBUG_FLAG 1

namespace Amanzi {

// -- Initialize owned (dependent) variables.
void MPCSubsurface::setup(const Teuchos::Ptr<State>& S) {
  // set up keys
  domain_name_ = plist_->get<std::string>("domain name", "domain");
  temp_key_ = getKey(domain_name_, "temperature");
  pres_key_ = getKey(domain_name_, "pressure");
  e_key_ = getKey(domain_name_, "energy");
  wc_key_ = getKey(domain_name_, "water_content");
  tc_key_ = getKey(domain_name_, "thermal_conductivity");
  uw_tc_key_ = getKey(domain_name_, "upwind_thermal_conductivity");
  kr_key_ = getKey(domain_name_, "relative_permeability");
  uw_kr_key_ = getKey(domain_name_, "upwind_relative_permeability");
  enth_key_ = getKey(domain_name_, "enthalpy");
  hkr_key_ = getKey(domain_name_, "enthalpy_times_relative_permeability");
  uw_hkr_key_ = getKey(domain_name_, "upwind_enthalpy_times_relative_permeability");
  energy_flux_key_ = getKey(domain_name_, "energy_flux");
  mass_flux_key_ = getKey(domain_name_, "mass_flux");
  mass_flux_dir_key_ = getKey(domain_name_, "mass_flux_direction");
  rho_key_ = getKey(domain_name_, "mass_density_liquid");

  // supress energy's vision of advective terms as we can do better
  Teuchos::Array<std::string> pk_order = plist_->get< Teuchos::Array<std::string> >("PKs order");
  if (plist_->isParameter("supress Jacobian terms: div hq / dp")) {
    Errors::Message message("MPC Incorrect input: parameter \"supress Jacobian terms: div hq / dp\" changed to \"supress Jacobian terms: div hq / dp,T\" for clarity.");
    Exceptions::amanzi_throw(message);
  }
  if (!plist_->get<bool>("supress Jacobian terms: div hq / dp,T", false)) {
    if (plist_->sublist("PKs").sublist(pk_order[1]).isParameter("supress advective terms in preconditioner")
        && !plist_->sublist("PKs").sublist(pk_order[1]).get("supress advective terms in preconditioner",false)) {
       Errors::Message msg("MPC Incorrect input: options \"supress Jacobian terms: div hq / dp,T\" and subsurface energy PK option \"supress advective terms in preconditioner\" should not both be false, as these include some of the same Jacobian information.\n Recommended: Enable/suppress the latter.");

       Exceptions::amanzi_throw(msg);
    }
  }
  
  // set up the sub-pks
  StrongMPC<PKPhysicalBDFBase>::setup(S);
  mesh_ = S->GetMesh(domain_name_);

  // set up debugger
  db_ = sub_pks_[0]->debugger();

  // Get the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Operator> pcA = sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Operator> pcB = sub_pks_[1]->preconditioner();
  if (pcA->DomainMap().HasComponent("face")) {
    is_fv_ = false;
  } else {
    is_fv_ = true;
  }
  
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
  } else if (precon_string == "picard") {
    precon_type_ = PRECON_PICARD;
  } else if (precon_string == "ewc") {
    ASSERT(0);
    precon_type_ = PRECON_EWC;
  } else if (precon_string == "smart ewc") {
    ASSERT(0);
    precon_type_ = PRECON_EWC;
  } else {
    Errors::Message message(std::string("Invalid preconditioner type ")+precon_string);
    Exceptions::amanzi_throw(message);
  }

  // create offdiagonal blocks
  if (precon_type_ != PRECON_NONE && precon_type_ != PRECON_BLOCK_DIAGONAL) {
    std::vector<AmanziMesh::Entity_kind> locations2(2);
    std::vector<std::string> names2(2);
    std::vector<int> num_dofs2(2,1);
    locations2[0] = AmanziMesh::CELL;
    names2[0] = "cell";

    // Create the block for derivatives of mass conservation with respect to temperature
    // -- derivatives of kr with respect to temperature
    if (!plist_->get<bool>("supress Jacobian terms: d div q / dT", false)) {
      // need to upwind dkr/dT
      if (!is_fv_) {
        Key dkrdT_key = getDerivKey(uw_kr_key_, temp_key_);
        S->RequireField(dkrdT_key, name_)
            ->SetMesh(mesh_)->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
        S->GetField(dkrdT_key,name_)->set_io_vis(false);
        upwinding_dkrdT_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                getDerivKey(kr_key_, temp_key_),
                dkrdT_key, mass_flux_dir_key_, 1.e-8));
      }

      // set up the operator
      Teuchos::ParameterList divq_plist(plist_->sublist("PKs").sublist(pk_order[0]).sublist("Diffusion PC"));
      if (is_fv_) divq_plist.set("Newton correction", "true Jacobian");
      else divq_plist.set("Newton correction", "approximate Jacobian");
      divq_plist.set("exclude primary terms", true);
      Operators::OperatorDiffusionFactory opfactory;
      ddivq_dT_ = opfactory.CreateWithGravity(divq_plist, mesh_);
      dWC_dT_block_ = ddivq_dT_->global_operator();
    }

    // -- derivatives of water content with respect to temperature
    if (dWC_dT_block_ == Teuchos::null) {
      dWC_dT_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, mesh_));
      dWC_dT_block_ = dWC_dT_->global_operator();
    } else {
      dWC_dT_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, dWC_dT_block_));
    }

    // Create the block for derivatives of energy conservation with respect to pressure
    // -- derivatives of thermal conductivity with respect to pressure
    if (!plist_->get<bool>("supress Jacobian terms: d div K grad T / dp", false)) {
      // need to upwind dKappa/dp
      if (!is_fv_) {
        Key uw_dKappa_dp_key = getDerivKey(uw_tc_key_, pres_key_);
        Key dKappa_dp_key = getDerivKey(tc_key_, pres_key_);
        S->RequireField(uw_dKappa_dp_key, name_)
            ->SetMesh(mesh_)->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
        S->GetField(uw_dKappa_dp_key,name_)->set_io_vis(false);

        upwinding_dKappa_dp_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                dKappa_dp_key, uw_dKappa_dp_key,
                energy_flux_key_, 1.e-8));
        // upwinding_dKappa_dp_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
        //         dKappa_dp_key, uw_dKappa_dp_key));
      }

      // set up the operator
      Teuchos::ParameterList ddivKgT_dp_plist(plist_->sublist("PKs").sublist(pk_order[1]).sublist("Diffusion PC"));
      if (is_fv_) ddivKgT_dp_plist.set("Newton correction", "true Jacobian");
      else ddivKgT_dp_plist.set("Newton correction", "approximate Jacobian");

      ddivKgT_dp_plist.set("exclude primary terms", true);
      Operators::OperatorDiffusionFactory opfactory;
      if (dE_dp_block_ == Teuchos::null) {
        ddivKgT_dp_ = opfactory.Create(ddivKgT_dp_plist, mesh_);
        dE_dp_block_ = ddivKgT_dp_->global_operator();
      } else {
        ddivKgT_dp_ = opfactory.Create(ddivKgT_dp_plist, dE_dp_block_);
      }
      ddivKgT_dp_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[0]->BCs());
    }


    // -- derivatives of advection term
    if (!plist_->get<bool>("supress Jacobian terms: div hq / dp,T", false)) {
      // derivative with respect to pressure
      Teuchos::ParameterList divhq_dp_plist(plist_->sublist("PKs").sublist(pk_order[0]).sublist("Diffusion PC"));
      if (is_fv_) divhq_dp_plist.set("Newton correction", "true Jacobian");
      else divhq_dp_plist.set("Newton correction", "approximate Jacobian");

      Operators::OperatorDiffusionFactory opfactory;
      if (dE_dp_block_ == Teuchos::null) {
        ddivhq_dp_ = opfactory.CreateWithGravity(divhq_dp_plist, mesh_);
        dE_dp_block_ = ddivhq_dp_->global_operator();
      } else {
        ddivhq_dp_ = opfactory.CreateWithGravity(divhq_dp_plist, dE_dp_block_);
      }

      // derivative with respect to temperature
      Teuchos::ParameterList divhq_dT_plist(plist_->sublist("PKs").sublist(pk_order[0]).sublist("Diffusion PC"));
      divhq_dT_plist.set("exclude primary terms", true);
      if (is_fv_) divhq_dT_plist.set("Newton correction", "true Jacobian");
      else divhq_dT_plist.set("Newton correction", "approximate Jacobian");

      ddivhq_dT_ = opfactory.CreateWithGravity(divhq_dT_plist, pcB);

      // need a field, evaluator, and upwinding for h * kr * rho/mu
      // -- first the evaluator
      Teuchos::ParameterList hkr_eval_list;
      hkr_eval_list.set("evaluator name", hkr_key_);
      Teuchos::Array<std::string> deps(2);
      deps[0] = enth_key_; deps[1] = kr_key_;
      hkr_eval_list.set("evaluator dependencies", deps);
      Teuchos::RCP<FieldEvaluator> hkr_eval =
          Teuchos::rcp(new Relations::MultiplicativeEvaluator(hkr_eval_list));

      // -- now the field
      names2[1] = "boundary_face";
      locations2[1] = AmanziMesh::BOUNDARY_FACE;
      S->RequireField(hkr_key_)->SetMesh(mesh_)->SetGhosted()
          ->AddComponents(names2, locations2, num_dofs2);
      S->SetFieldEvaluator(hkr_key_, hkr_eval);

      S->RequireField(uw_hkr_key_, name_)
          ->SetMesh(mesh_)->SetGhosted()
          ->SetComponent("face", AmanziMesh::FACE, 1);
      S->GetField(uw_hkr_key_,name_)->set_io_vis(false);

      std::string method_name = plist_->sublist("PKs").sublist(pk_order[0])
          .get<std::string>("relative permeability method", "upwind with gravity");
      if (method_name != "upwind with Darcy flux") {
        Errors::Message msg;
        msg << "Subsurface coupler with advective Jacobian terms only supports a Richards upwind scheme of "
            << "\"upwind with Darcy flux\", but the method \"" << method_name << "\" was requested.";
        Exceptions::amanzi_throw(msg);
      }
      upwinding_hkr_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
              hkr_key_, uw_hkr_key_, mass_flux_dir_key_, 1.e-8));
      
      if (!is_fv_) {
        // -- and the upwinded field
        S->RequireField(getDerivKey(uw_hkr_key_, pres_key_), name_)
            ->SetMesh(mesh_)->SetGhosted()
            ->SetComponent("face", AmanziMesh::FACE, 1);
        S->GetField(getDerivKey(uw_hkr_key_, pres_key_),name_)
            ->set_io_vis(false);
        S->RequireField(getDerivKey(uw_hkr_key_, temp_key_), name_)
            ->SetMesh(mesh_)->SetGhosted()
            ->SetComponent("face", AmanziMesh::FACE, 1);
        S->GetField(getDerivKey(uw_hkr_key_, temp_key_), name_)
            ->set_io_vis(false);

        // -- and the upwinding
        upwinding_dhkr_dp_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                getDerivKey(hkr_key_, pres_key_),
                getDerivKey(uw_hkr_key_, pres_key_),
                mass_flux_dir_key_, 1.e-8));
        upwinding_dhkr_dT_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                getDerivKey(hkr_key_, temp_key_),
                getDerivKey(uw_hkr_key_, temp_key_),
                mass_flux_dir_key_, 1.e-8));
      }
    }

    // -- derivatives of energy with respect to pressure
    if (dE_dp_block_ == Teuchos::null) {
      dE_dp_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, mesh_));
      dE_dp_block_ = dE_dp_->global_operator();
    } else {
      dE_dp_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, dE_dp_block_));
    }

    ASSERT(dWC_dT_block_ != Teuchos::null);
    ASSERT(dE_dp_block_ != Teuchos::null);
    preconditioner_->SetOperatorBlock(0, 1, dWC_dT_block_);
    preconditioner_->SetOperatorBlock(1, 0, dE_dp_block_);
  }
  
  // set up sparsity structure
  preconditioner_->SymbolicAssembleMatrix();

  // create the linear solver
  if (plist_->isSublist("linear solver")) {
    Teuchos::ParameterList& lin_solver_list = plist_->sublist("linear solver");
    AmanziSolvers::LinearOperatorFactory<Operators::TreeOperator, TreeVector, TreeVectorSpace> fac;
    linsolve_preconditioner_ = fac.Create(lin_solver_list, preconditioner_);
  } else {
    linsolve_preconditioner_ = preconditioner_;
  }
  
  // create the EWC delegate
  if (plist_->isSublist("ewc delegate")) {
    Teuchos::RCP<Teuchos::ParameterList> sub_ewc_list = Teuchos::sublist(plist_, "ewc delegate");
    sub_ewc_list->set("PK name", name_);
    sub_ewc_list->set("domain key", domain_name_);
    ewc_ = Teuchos::rcp(new MPCDelegateEWCSubsurface(*sub_ewc_list));
    Teuchos::RCP<PermafrostModel> model = Teuchos::rcp(new PermafrostModel());
    ewc_->set_model(model);
    ewc_->setup(S);
  }    
}

void MPCSubsurface::initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC<PKPhysicalBDFBase>::initialize(S);
  if (ewc_ != Teuchos::null) ewc_->initialize(S);

  // initialize offdiagonal operators
  richards_pk_ = Teuchos::rcp_dynamic_cast<Flow::Richards>(sub_pks_[0]);
  ASSERT(richards_pk_ != Teuchos::null);

  if (ddivq_dT_ != Teuchos::null) {
    if (!is_fv_) {
      Key dkrdT_key = getDerivKey(uw_kr_key_, temp_key_);
      S->GetFieldData(dkrdT_key,name_)->PutScalar(1.0);
      S->GetField(dkrdT_key,name_)->set_initialized();
    }

    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    AmanziGeometry::Point g(3);
    g[0] = (*gvec)[0]; g[1] = (*gvec)[1]; g[2] = (*gvec)[2];
    ddivq_dT_->SetGravity(g);    
    ddivq_dT_->SetBCs(sub_pks_[0]->BCs(), sub_pks_[1]->BCs());
    ddivq_dT_->SetTensorCoefficient(richards_pk_->K_);
  }

  if (ddivKgT_dp_ != Teuchos::null) {
    if (!is_fv_) {
      Key uw_dKappa_dp_key = getDerivKey(uw_tc_key_, pres_key_);
      S->GetFieldData(uw_dKappa_dp_key,name_)->PutScalar(1.0);
      S->GetField(uw_dKappa_dp_key,name_)->set_initialized();
    }

    ddivKgT_dp_->SetTensorCoefficient(Teuchos::null);
  }

  if (ddivhq_dp_ != Teuchos::null) {

    S->GetFieldData(uw_hkr_key_, name_)->PutScalar(1.);
    S->GetField(uw_hkr_key_, name_)->set_initialized();
    if (!is_fv_) {
      S->GetFieldData(getDerivKey(uw_hkr_key_, pres_key_), name_)->PutScalar(1.);
      S->GetField(getDerivKey(uw_hkr_key_, pres_key_), name_)->set_initialized();
      S->GetFieldData(getDerivKey(uw_hkr_key_, temp_key_), name_)->PutScalar(1.);
      S->GetField(getDerivKey(uw_hkr_key_, temp_key_), name_)->set_initialized();
    }

    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    AmanziGeometry::Point g(3);
    g[0] = (*gvec)[0]; g[1] = (*gvec)[1]; g[2] = (*gvec)[2];
    ddivhq_dp_->SetGravity(g);    
    ddivhq_dp_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[0]->BCs());
    ddivhq_dp_->SetTensorCoefficient(richards_pk_->K_);

    ddivhq_dT_->SetGravity(g);    
    ddivhq_dT_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[1]->BCs());
    ddivhq_dT_->SetTensorCoefficient(richards_pk_->K_);
  }  

}


void MPCSubsurface::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  StrongMPC<PKPhysicalBDFBase>::set_states(S,S_inter,S_next);
  if (ewc_ != Teuchos::null) ewc_->set_states(S,S_inter,S_next);
}

void MPCSubsurface::commit_state(double dt, const Teuchos::RCP<State>& S) {
  StrongMPC<PKPhysicalBDFBase>::commit_state(dt,S);
  if (ewc_ != Teuchos::null) ewc_->commit_state(dt,S);
}

// update the predictor to be physically consistent
bool MPCSubsurface::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> up0,
        Teuchos::RCP<TreeVector> up) {
  bool modified(false);
  if (ewc_ != Teuchos::null) {
    modified = ewc_->ModifyPredictor(h, up);
    if (modified) ChangedSolution();
  }

  // potentially update faces
  modified |= StrongMPC<PKPhysicalBDFBase>::ModifyPredictor(h, up0, up);
  return modified;
}


// updates the preconditioner
void MPCSubsurface::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h, bool assemble) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (precon_type_ == PRECON_NONE) {
    // nothing to do
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    StrongMPC<PKPhysicalBDFBase>::UpdatePreconditioner(t,up,h);
  } else if (precon_type_ == PRECON_PICARD || precon_type_ == PRECON_EWC) {
    StrongMPC<PKPhysicalBDFBase>::UpdatePreconditioner(t,up,h);

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
      if (is_fv_) {
        dkrdT = S_next_->GetFieldData(getDerivKey(kr_key_, temp_key_));
      } else {
        S_next_->GetFieldData(getDerivKey(uw_kr_key_, temp_key_), name_)
            ->PutScalar(0.);
        upwinding_dkrdT_->Update(S_next_.ptr());
        dkrdT = S_next_->GetFieldData(getDerivKey(uw_kr_key_, temp_key_));
      }

      // form the operator
      Teuchos::RCP<const CompositeVector> kr_uw = S_next_->GetFieldData(uw_kr_key_);
      Teuchos::RCP<const CompositeVector> flux = S_next_->GetFieldData(mass_flux_key_);
      Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(rho_key_);

      ddivq_dT_->SetDensity(rho);
      ddivq_dT_->SetScalarCoefficient(kr_uw, dkrdT);
      ddivq_dT_->UpdateMatrices(flux.ptr(),
              up->SubVector(0)->Data().ptr());              
      ddivq_dT_->UpdateMatricesNewtonCorrection(flux.ptr(),
              up->SubVector(0)->Data().ptr());              
      ddivq_dT_->ApplyBCs(false, true);
    }

    // -- dWC/dT diagonal term
    S_next_->GetFieldEvaluator(wc_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, temp_key_);
    Teuchos::RCP<const CompositeVector> dWC_dT =
      S_next_->GetFieldData(getDerivKey(wc_key_, temp_key_));
    dWC_dT_->AddAccumulationTerm(*dWC_dT->ViewComponent("cell", false), h);

    // std::cout << "1/h * DWC/DT" << std::endl;
    // CompositeVector dbg(*dWC_dT); dbg = *dWC_dT; dbg.Scale(1./h); dbg.Print(std::cout);
    

    // dE / dp block
    // -- d Kappa / dp
    if (ddivKgT_dp_ != Teuchos::null) {
      // Update and upwind thermal conductivity
      S_next_->GetFieldEvaluator(tc_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, pres_key_);

      Teuchos::RCP<const CompositeVector> dKappa_dp;
      if (is_fv_) {
        dKappa_dp = S_next_->GetFieldData(getDerivKey(tc_key_, pres_key_));
      } else {
        S_next_->GetFieldData(getDerivKey(uw_tc_key_, pres_key_), name_)
            ->PutScalar(0.);
        upwinding_dKappa_dp_->Update(S_next_.ptr(), db_.ptr());
        dKappa_dp = S_next_->GetFieldData(getDerivKey(uw_tc_key_, pres_key_));
      }

      // form the operator
      Teuchos::RCP<const CompositeVector> uw_Kappa =
        S_next_->GetFieldData(uw_tc_key_);
      Teuchos::RCP<const CompositeVector> flux =
        S_next_->GetFieldData(energy_flux_key_);
      ddivKgT_dp_->SetScalarCoefficient(uw_Kappa, dKappa_dp);
      ddivKgT_dp_->UpdateMatrices(flux.ptr(),
              up->SubVector(1)->Data().ptr());
      ddivKgT_dp_->UpdateMatricesNewtonCorrection(flux.ptr(),
              up->SubVector(1)->Data().ptr());
      ddivKgT_dp_->ApplyBCs(false, true);
    }

    // -- d adv / dp   This one is a bit more complicated...
    if (ddivhq_dp_ != Teuchos::null) {
      // Update and upwind enthalpy * kr * rho/mu
      // -- update values
      S_next_->GetFieldEvaluator(hkr_key_)
          ->HasFieldChanged(S_next_.ptr(), name_);
      S_next_->GetFieldEvaluator(hkr_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, pres_key_);
      S_next_->GetFieldEvaluator(hkr_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, temp_key_);

      Teuchos::RCP<const CompositeVector> denth_kr_dp_uw;
      Teuchos::RCP<const CompositeVector> denth_kr_dT_uw;

      Teuchos::RCP<const CompositeVector> enth_kr =
        S_next_->GetFieldData(hkr_key_);
      Teuchos::RCP<CompositeVector> enth_kr_uw =
        S_next_->GetFieldData(uw_hkr_key_, name_);
      enth_kr_uw->PutScalar(0.);

      enth_kr_uw->ViewComponent("face",false)
          ->Export(*enth_kr->ViewComponent("boundary_face",false),
                   mesh_->exterior_face_importer(), Insert);
      upwinding_hkr_->Update(S_next_.ptr(), db_.ptr());

      if (richards_pk_->clobber_surf_kr_) {
        // -- stick zeros in the boundary faces
        Epetra_MultiVector enth_kr_bf(*enth_kr->ViewComponent("boundary_face",false));
        enth_kr_bf.PutScalar(0.0);
        enth_kr_uw->ViewComponent("face",false)->Export(enth_kr_bf,
                mesh_->exterior_face_importer(), Insert);
      }
      
      if (is_fv_) {
        denth_kr_dp_uw = S_next_->GetFieldData(getDerivKey(hkr_key_, pres_key_));
        denth_kr_dT_uw = S_next_->GetFieldData(getDerivKey(hkr_key_, temp_key_));
      } else {

        Teuchos::RCP<const CompositeVector> denth_kr_dp =
            S_next_->GetFieldData(getDerivKey(hkr_key_, pres_key_));
        Teuchos::RCP<const CompositeVector> denth_kr_dT =
            S_next_->GetFieldData(getDerivKey(hkr_key_, temp_key_));

        // -- zero target data (may be unnecessary?)
        Teuchos::RCP<CompositeVector> denth_kr_dp_uw_nc =
            S_next_->GetFieldData(getDerivKey(uw_hkr_key_, pres_key_), name_);
        denth_kr_dp_uw_nc->PutScalar(0.);
        Teuchos::RCP<CompositeVector> denth_kr_dT_uw_nc =
            S_next_->GetFieldData(getDerivKey(uw_hkr_key_, temp_key_), name_);
        denth_kr_dT_uw_nc->PutScalar(0.);

        // -- copy boundary faces into upwinded vector
        denth_kr_dp_uw_nc->ViewComponent("face",false)
            ->Export(*denth_kr_dp->ViewComponent("boundary_face",false),
                     mesh_->exterior_face_importer(), Insert);
        denth_kr_dT_uw_nc->ViewComponent("face",false)
            ->Export(*denth_kr_dT->ViewComponent("boundary_face",false),
                     mesh_->exterior_face_importer(), Insert);

        // -- upwind      
        upwinding_dhkr_dp_->Update(S_next_.ptr(), db_.ptr());
        upwinding_dhkr_dT_->Update(S_next_.ptr(), db_.ptr());

        // -- clobber
        if (richards_pk_->clobber_surf_kr_) {
          // -- stick zeros in the boundary faces
          Epetra_MultiVector enth_kr_bf(*enth_kr->ViewComponent("boundary_face",false));
          enth_kr_bf.PutScalar(0.0);
          denth_kr_dp_uw_nc->ViewComponent("face",false)->Export(enth_kr_bf,
                  mesh_->exterior_face_importer(), Insert);
          denth_kr_dT_uw_nc->ViewComponent("face",false)->Export(enth_kr_bf,
                  mesh_->exterior_face_importer(), Insert);
        }

        denth_kr_dp_uw =
            S_next_->GetFieldData(getDerivKey(uw_hkr_key_, pres_key_));
        denth_kr_dT_uw =
            S_next_->GetFieldData(getDerivKey(uw_hkr_key_, temp_key_));

      }
      
      Teuchos::RCP<const CompositeVector> flux = S_next_->GetFieldData(mass_flux_key_);
      Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(rho_key_);

      // form the operator: pressure component
      ddivhq_dp_->SetDensity(rho);
      ddivhq_dp_->SetScalarCoefficient(enth_kr_uw, denth_kr_dp_uw);
      // -- update the local matrices, div h * kr grad
      ddivhq_dp_->UpdateMatrices(Teuchos::null, Teuchos::null);
      // -- determine the advective fluxes, q_a = h * kr grad p
      CompositeVector adv_flux(*flux, INIT_MODE_ZERO);
      Teuchos::Ptr<const CompositeVector> adv_flux_ptr(&adv_flux);
      ddivhq_dp_->UpdateFlux(*up->SubVector(0)->Data(), adv_flux);
      // -- add in components div (d h*kr / dp) grad q_a / (h*kr)
      ddivhq_dp_->UpdateMatricesNewtonCorrection(adv_flux_ptr, up->SubVector(0)->Data().ptr());
      ddivhq_dp_->ApplyBCs(false, true);

      // form the operator: temperature component
      ddivhq_dT_->SetDensity(rho);
      ddivhq_dT_->SetScalarCoefficient(enth_kr_uw, denth_kr_dT_uw);
      // -- add in components div (d h*kr / dp) grad q_a / (h*kr)
      ddivhq_dT_->UpdateMatrices(adv_flux_ptr, up->SubVector(0)->Data().ptr());
      ddivhq_dT_->UpdateMatricesNewtonCorrection(adv_flux_ptr, up->SubVector(0)->Data().ptr());
      ddivhq_dT_->ApplyBCs(false, true);

    }

    // -- dE/dp diagonal term
    S_next_->GetFieldEvaluator(e_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, pres_key_);
    Teuchos::RCP<const CompositeVector> dE_dp =
      S_next_->GetFieldData(getDerivKey(e_key_, pres_key_));
    dE_dp_->AddAccumulationTerm(*dE_dp->ViewComponent("cell", false), h);

    // std::cout << "1/h * DE/Dp" << std::endl;
    // dbg = *dE_dp; dbg.Scale(1./h); dbg.Print(std::cout);
    
    // write for debugging
    std::vector<std::string> vnames;
    vnames.push_back("  dwc_dT"); vnames.push_back("  de_dp"); 
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(dWC_dT.ptr()); vecs.push_back(dE_dp.ptr());
    db_->WriteVectors(vnames, vecs, false);

    // finally assemble the full system, dump if requested, and form the inverse
    if (assemble) {
      preconditioner_->AssembleMatrix();
      if (dump_) {
        std::stringstream filename;
        filename << "Subsurface_PC_" << S_next_->cycle() << ".txt";
        EpetraExt::RowMatrixToMatlabFile(filename.str().c_str(), *preconditioner_->A());
      }
      Teuchos::ParameterList& pc_sublist = plist_->sublist("preconditioner");
      preconditioner_->InitPreconditioner(pc_sublist);
    }
  }
  
  if (precon_type_ == PRECON_EWC) {
    ewc_->UpdatePreconditioner(t,up,h);
  }
}


// -----------------------------------------------------------------------------
// Wrapper to call the requested preconditioner.
// -----------------------------------------------------------------------------
int MPCSubsurface::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
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
  
  int ierr;
  if (precon_type_ == PRECON_NONE) {
    *Pu = *u;
    ierr = 1;
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    ierr = StrongMPC<PKPhysicalBDFBase>::ApplyPreconditioner(u,Pu);
  } else if (precon_type_ == PRECON_PICARD) {
    ierr = linsolve_preconditioner_->ApplyInverse(*u, *Pu);
  } else if (precon_type_ == PRECON_EWC) {
    ierr = linsolve_preconditioner_->ApplyInverse(*u, *Pu);

  //   if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //     *vo_->os() << "PC_std * residuals:" << std::endl;
  //     std::vector<std::string> vnames;
  //     vnames.push_back("  PC*r_p"); vnames.push_back("  PC*r_T"); 
  //     std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  //     vecs.push_back(Pu->SubVector(0)->Data().ptr()); 
  //     vecs.push_back(Pu->SubVector(1)->Data().ptr()); 
  //     db_->WriteVectors(vnames, vecs, true);
  //   }

  //   // make sure we can back-calc face corrections that preserve residuals on faces
  //   Teuchos::RCP<TreeVector> res0 = Teuchos::rcp(new TreeVector(*u));
  //   res0->PutScalar(0.);
  //   Teuchos::RCP<TreeVector> Pu_std = Teuchos::rcp(new TreeVector(*Pu));
  //   *Pu_std = *Pu;

  //   // call EWC, which does Pu_p <-- Pu_p_std + dPu_p
  //   ewc_->ApplyPreconditioner(u, Pu);

  //   // calculate dPu_lambda from dPu_p
  //   Pu_std->Update(1.0, *Pu, -1.0);
  //   mfd_preconditioner_->UpdateConsistentFaceCorrection(*res0, Pu_std.ptr());

  //   // update Pu_lambda <-- Pu_lambda_std + dPu_lambda
  //   Pu->SubVector(0)->Data()->ViewComponent("face",false)->Update(1.,
  //           *Pu_std->SubVector(0)->Data()->ViewComponent("face",false), 1.);
  //   Pu->SubVector(1)->Data()->ViewComponent("face",false)->Update(1.,
  //           *Pu_std->SubVector(1)->Data()->ViewComponent("face",false), 1.);

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


// AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
//     MPCSubsurface::ModifyCorrection(double h,
//                                     Teuchos::RCP<const TreeVector> res,
//                                     Teuchos::RCP<const TreeVector> u,
//                                     Teuchos::RCP<TreeVector> du) {

//   if (vo_->os_OK(Teuchos::VERB_HIGH)) {
//     *vo_->os() << "NKA * PC * residuals:" << std::endl;
//     std::vector<std::string> vnames;
//     vnames.push_back("  NKA*PC*r_p"); vnames.push_back("  NKA*PC*r_T"); 
//     std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
//     vecs.push_back(du->SubVector(0)->Data().ptr()); 
//     vecs.push_back(du->SubVector(1)->Data().ptr()); 
//     db_->WriteVectors(vnames, vecs, true);
//   }

//   // if (precon_type_ == PRECON_EWC) {
//   //   // make sure we can back-calc face corrections that preserve residuals on faces
//   //   Teuchos::RCP<TreeVector> res0 = Teuchos::rcp(new TreeVector(*res));
//   //   res0->PutScalar(0.);
//   //   Teuchos::RCP<TreeVector> du_std = Teuchos::rcp(new TreeVector(*du));
//   //   *du_std = *du;

//   //   // call EWC, which does du_p <-- du_p_std + ddu_p
//   //   ewc_->ApplyPreconditioner(res, du);

//   //   // calculate ddu_lambda from ddu_p
//   //   du_std->Update(1.0, *du, -1.0);
//   //   preconditioner_->UpdateConsistentFaceCorrection(*res0, du_std.ptr());

//   //   // update du_lambda <-- du_lambda_std + ddu_lambda
//   //   du->SubVector(0)->Data()->ViewComponent("face",false)->Update(1.,
//   //           *du_std->SubVector(0)->Data()->ViewComponent("face",false), 1.);
//   //   du->SubVector(1)->Data()->ViewComponent("face",false)->Update(1.,
//   //           *du_std->SubVector(1)->Data()->ViewComponent("face",false), 1.);
//   // }

//   return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
// }

} // namespace
