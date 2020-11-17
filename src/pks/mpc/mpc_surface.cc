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

#include "surface_ice_model.hh"

#include "mpc_delegate_ewc_surface.hh"
#include "mpc_surface.hh"

#define DEBUG_FLAG 1

namespace Amanzi {

// -- Initialize owned (dependent) variables.
void MPCSurface::Setup(const Teuchos::Ptr<State>& S)
{
  auto pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
  domain_ = plist_->get<std::string>("domain name");

  temp_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
  pres_key_ = Keys::readKey(*plist_, domain_, "pressure", "pressure");
  e_key_ = Keys::readKey(*plist_, domain_, "energy", "energy");
  wc_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");

  kr_key_ = Keys::readKey(*plist_, domain_, "overland conductivity", "overland_conductivity");
  kr_uw_key_ = Keys::readKey(*plist_, domain_, "upwind overland conductivity", "upwind_overland_conductivity");
  potential_key_ = Keys::readKey(*plist_, domain_, "potential", "pres_elev");
  pd_bar_key_ = Keys::readKey(*plist_, domain_, "ponded depth, negative", "ponded_depth_bar");
  mass_flux_key_ = Keys::readKey(*plist_, domain_, "mass flux", "mass_flux");

  // make sure the overland flow pk does not rescale the preconditioner -- we want it in h
  pks_list_->sublist(pk_order[0]).set("scale preconditioner to pressure", false);

  // set up the sub-pks
  StrongMPC<PK_PhysicalBDF_Default>::Setup(S);
  mesh_ = S->GetMesh(domain_);

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
  preconditioner_->set_operator_block(0, 0, pcA);
  preconditioner_->set_operator_block(1, 1, pcB);

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
      divq_plist.set("include Newton correction", true);
      divq_plist.set("exclude primary terms", true);
      Operators::PDE_DiffusionFactory opfactory;
      ddivq_dT_ = opfactory.Create(divq_plist, mesh_);
      dWC_dT_block_ = ddivq_dT_->global_operator();
    }

    // -- derivatives of water content with respect to temperature are zero on
    // -- the surface.

    // -- derivatives of energy with respect to pressure
    if (dE_dp_block_ == Teuchos::null) {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
      dE_dp_block_ = dE_dp_->global_operator();
    } else {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, dE_dp_block_));
    }

    preconditioner_->set_operator_block(0, 1, dWC_dT_block_);
    preconditioner_->set_operator_block(1, 0, dE_dp_block_);
    preconditioner_->set_inverse_parameters(plist_->sublist("inverse"));
  }

  // create the EWC delegate
  if (plist_->isSublist("surface ewc delegate")) {
    Teuchos::RCP<Teuchos::ParameterList> surf_ewc_list = Teuchos::sublist(plist_, "surface ewc delegate");
    surf_ewc_list->set("PK name", name_);
    surf_ewc_list->set("domain name", domain_);
    ewc_ = Teuchos::rcp(new MPCDelegateEWCSurface(*surf_ewc_list));

    Teuchos::RCP<EWCModelBase> model = Teuchos::rcp(new SurfaceIceModel());
    ewc_->set_model(model);
    ewc_->setup(S);
  }
}


void MPCSurface::Initialize(const Teuchos::Ptr<State>& S) {
  StrongMPC<PK_PhysicalBDF_Default>::Initialize(S);
  if (ewc_ != Teuchos::null) ewc_->initialize(S);

  if (ddivq_dT_ != Teuchos::null) {
    ddivq_dT_->SetBCs(sub_pks_[0]->BCs(), sub_pks_[1]->BCs());
    ddivq_dT_->SetTensorCoefficient(Teuchos::null);
  }
}


void MPCSurface::set_states(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next)
{
  StrongMPC<PK_PhysicalBDF_Default>::set_states(S,S_inter,S_next);
  if (ewc_ != Teuchos::null) ewc_->set_states(S,S_inter,S_next);
}


void MPCSurface::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  StrongMPC<PK_PhysicalBDF_Default>::CommitStep(t_old, t_new, S);
  if (ewc_ != Teuchos::null) {
    double dt = t_new - t_old;
    ewc_->commit_state(dt,S);
  }
}


// update the predictor to be physically consistent
bool MPCSurface::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> up0,
        Teuchos::RCP<TreeVector> up)
{
  bool modified(false);
  if (ewc_ != Teuchos::null) {
    modified = ewc_->ModifyPredictor(h, up);
    if (modified) ChangedSolution();
  }

  // potentially update faces
  modified |= StrongMPC<PK_PhysicalBDF_Default>::ModifyPredictor(h, up0, up);
  return modified;
}


// updates the preconditioner
void MPCSurface::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (precon_type_ == PRECON_NONE) {
    // nothing to do
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    StrongMPC::UpdatePreconditioner(t,up,h);
  } else if (precon_type_ == PRECON_PICARD || precon_type_ == PRECON_EWC) {
    preconditioner_->InitOffdiagonals(); // zero out offdiagonal blocks
    StrongMPC::UpdatePreconditioner(t,up,h);

    // dWC / dT block
    // -- dkr/dT
    if (ddivq_dT_ != Teuchos::null) {
      // -- update and upwind d kr / dT
      S_next_->GetFieldEvaluator(kr_key_)
          ->HasFieldDerivativeChanged(S_next_.ptr(), name_, temp_key_);
      Teuchos::RCP<const CompositeVector> dkrdT =
        S_next_->GetFieldData(Keys::getDerivKey(kr_key_, temp_key_));
      Teuchos::RCP<const CompositeVector> kr_uw =
        S_next_->GetFieldData(kr_uw_key_);
      Teuchos::RCP<const CompositeVector> flux =
        S_next_->GetFieldData(mass_flux_key_);

      S_next_->GetFieldEvaluator(potential_key_)
        ->HasFieldChanged(S_next_.ptr(), name_);
      Teuchos::RCP<const CompositeVector> pres_elev =
        S_next_->GetFieldData(potential_key_);

      // form the operator
      ddivq_dT_->SetScalarCoefficient(kr_uw, dkrdT);
      ddivq_dT_->UpdateMatrices(flux.ptr(), pres_elev.ptr());
      ddivq_dT_->UpdateMatricesNewtonCorrection(flux.ptr(), pres_elev.ptr());
      ddivq_dT_->ApplyBCs(false, true, false);
    }

    // -- dE/dp diagonal term
    S_next_->GetFieldEvaluator(e_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, pres_key_);
    Teuchos::RCP<const CompositeVector> dE_dp =
      S_next_->GetFieldData(Keys::getDerivKey(e_key_, pres_key_));

    // -- scale to dE/dh
    S_next_->GetFieldEvaluator(pd_bar_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, pres_key_);
    auto dh_dp = S_next_->GetFieldData(Keys::getDerivKey(pd_bar_key_, pres_key_));

    // -- add it in
    CompositeVector dE_dh(dE_dp->Map());
    dE_dh.ReciprocalMultiply(1./h, *dh_dp, *dE_dp, 0.);
    db_->WriteVector("  de_dp", dE_dp.ptr(), false);
    dE_dp_->AddAccumulationTerm(dE_dh, "cell");

    // write for debugging
    db_->WriteVector("  de_dp", dE_dp.ptr(), false);
    db_->WriteVector("  de_dh", Teuchos::ptr(&dE_dh), false);
  }
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
  } else if (precon_type_ == PRECON_PICARD || (precon_type_ == PRECON_EWC)) {
    ierr = preconditioner_->ApplyInverse(*u, *Pu);

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "PC * residuals:" << std::endl;
      std::vector<std::string> vnames;
      vnames.push_back("  PC*r_h"); vnames.push_back("  PC*r_T");
      std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
      vecs.push_back(Pu->SubVector(0)->Data().ptr());
      vecs.push_back(Pu->SubVector(1)->Data().ptr());
      db_->WriteVectors(vnames, vecs, true);
    }

    // tack on the variable change from h to p
    const Epetra_MultiVector& dh_dp =
      *S_next_->GetFieldData(Keys::getDerivKey(pd_bar_key_,pres_key_))->ViewComponent("cell",false);
    Epetra_MultiVector& Pu_c = *Pu->SubVector(0)->Data()->ViewComponent("cell",false);

    for (unsigned int c=0; c!=Pu_c.MyLength(); ++c) {
      Pu_c[0][c] /= dh_dp[0][c];
    }
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
