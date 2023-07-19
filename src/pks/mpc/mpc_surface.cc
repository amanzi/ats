/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

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

MPCSurface::MPCSurface(Teuchos::ParameterList& pk_tree_list,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree_list, global_list, S, soln),
    StrongMPC<PK_PhysicalBDF_Default>(pk_tree_list, global_list, S, soln)
{
  auto pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
  domain_ = plist_->get<std::string>("domain name");

  temp_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
  pres_key_ = Keys::readKey(*plist_, domain_, "pressure", "pressure");
  e_key_ = Keys::readKey(*plist_, domain_, "energy", "energy");
  wc_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");

  kr_key_ = Keys::readKey(*plist_, domain_, "overland conductivity", "overland_conductivity");
  kr_uw_key_ =
    Keys::readKey(*plist_, domain_, "upwind overland conductivity", "upwind_overland_conductivity");
  potential_key_ = Keys::readKey(*plist_, domain_, "potential", "pres_elev");
  pd_bar_key_ = Keys::readKey(*plist_, domain_, "ponded depth, negative", "ponded_depth_bar");
  water_flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");

  dump_ = plist_->get<bool>("dump preconditioner", false);

  // make sure the overland flow pk does not rescale the preconditioner -- we want it in h
  pks_list_->sublist(pk_order[0]).set("scale preconditioner to pressure", false);
}

// -- Initialize owned (dependent) variables.
void
MPCSurface::Setup()
{
  // set up the sub-pks
  StrongMPC<PK_PhysicalBDF_Default>::Setup();
  mesh_ = S_->GetMesh(domain_);

  // set up debugger
  db_ = sub_pks_[0]->debugger();

  // require these in case the PK did not do so already
  S_->Require<CompositeVector, CompositeVectorSpace>(pd_bar_key_, tag_next_)
    .SetMesh(S_->GetMesh(domain_))
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(pd_bar_key_, tag_next_);

  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    pd_bar_key_, tag_next_, pres_key_, tag_next_);

  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    e_key_, tag_next_, pres_key_, tag_next_);

  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    wc_key_, tag_next_, temp_key_, tag_next_);

  // Get the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Operator> pcA = sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Operator> pcB = sub_pks_[1]->preconditioner();

  // Create the combined operator
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcA->DomainMap()))));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcB->DomainMap()))));

  preconditioner_ =
    Teuchos::rcp(new Operators::TreeOperator(tvs, plist_->sublist("operator preconditioner")));
  preconditioner_->set_operator_block(0, 0, pcA);
  preconditioner_->set_operator_block(1, 1, pcB);

  // select the method used for preconditioning
  std::string precon_string = plist_->get<std::string>("preconditioner type", "picard");
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
    Errors::Message message(std::string("Invalid preconditioner type ") + precon_string);
    Exceptions::amanzi_throw(message);
  }

  // create offdiagonal blocks
  if (precon_type_ != PRECON_NONE && precon_type_ != PRECON_BLOCK_DIAGONAL) {
    // Create the block for derivatives of mass conservation with respect to temperature
    // -- derivatives of kr with respect to temperature
    if (precon_type_ != PRECON_NO_FLOW_COUPLING &&
        !plist_->get<bool>("supress Jacobian terms: d div surface q / dT", false)) {
      // set up the operator
      auto pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
      Teuchos::ParameterList divq_plist(
        pks_list_->sublist(pk_order[0]).sublist("diffusion preconditioner"));
      divq_plist.set("include Newton correction", true);
      divq_plist.set("exclude primary terms", true);
      Operators::PDE_DiffusionFactory opfactory;
      ddivq_dT_ = opfactory.Create(divq_plist, mesh_);
      dWC_dT_block_ = ddivq_dT_->global_operator();

      // require the derivative
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        kr_key_, tag_next_, temp_key_, tag_next_);
    }

    // -- derivatives of water content with respect to temperature are zero on
    // -- the surface.

    // -- derivatives of energy with respect to pressure
    if (dE_dp_block_ == Teuchos::null) {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
      dE_dp_block_ = dE_dp_->global_operator();
    } else {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, dE_dp_block_));
    }

    preconditioner_->set_operator_block(0, 1, dWC_dT_block_);
    preconditioner_->set_operator_block(1, 0, dE_dp_block_);
    preconditioner_->set_inverse_parameters(plist_->sublist("inverse"));
  }

  // This is currently broken, and shouldn't be used outside of this MPC... See #122
  // create the EWC delegate
  if (plist_->isSublist("surface ewc delegate")) {
    Teuchos::RCP<Teuchos::ParameterList> surf_ewc_list =
      Teuchos::sublist(plist_, "surface ewc delegate");
    surf_ewc_list->set("PK name", name_);
    surf_ewc_list->set("domain name", domain_);
    ewc_ = Teuchos::rcp(new MPCDelegateEWCSurface(*surf_ewc_list, S_));
    ewc_->set_tags(tag_current_, tag_next_);
    Teuchos::RCP<EWCModelBase> model = Teuchos::rcp(new SurfaceIceModel());
    ewc_->set_model(model);
    ewc_->setup();
  }
}


void
MPCSurface::Initialize()
{
  StrongMPC<PK_PhysicalBDF_Default>::Initialize();
  if (ewc_ != Teuchos::null) ewc_->initialize();

  if (ddivq_dT_ != Teuchos::null) {
    ddivq_dT_->SetBCs(sub_pks_[0]->BCs(), sub_pks_[1]->BCs());
    ddivq_dT_->SetTensorCoefficient(Teuchos::null);
  }
}

//void MPCSurface::set_tags(const Tag& tag_current, const Tag& tag_next)
//{
//  StrongMPC<PK_PhysicalBDF_Default>::set_tags(tag_current, tag_next);
//  if (ewc_ != Teuchos::null) ewc_->set_tags(tag_current, tag_next);
//}

void
MPCSurface::CommitStep(double t_old, double t_new, const Tag& tag)
{
  if (ewc_ != Teuchos::null) { ewc_->commit_state(); }
  StrongMPC<PK_PhysicalBDF_Default>::CommitStep(t_old, t_new, tag);
}


// update the predictor to be physically consistent
bool
MPCSurface::ModifyPredictor(double h,
                            Teuchos::RCP<const TreeVector> up0,
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
void
MPCSurface::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  if (precon_type_ == PRECON_NONE) {
    // nothing to do
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    StrongMPC::UpdatePreconditioner(t, up, h);
  } else if (precon_type_ == PRECON_PICARD || precon_type_ == PRECON_EWC) {
    preconditioner_->InitOffdiagonals(); // zero out offdiagonal blocks
    StrongMPC::UpdatePreconditioner(t, up, h);

    // dWC / dT block
    // -- dkr/dT
    if (ddivq_dT_ != Teuchos::null) {
      // -- update and upwind d kr / dT
      S_->GetEvaluator(kr_key_, tag_next_).UpdateDerivative(*S_, name_, temp_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> dkrdT =
        S_->GetDerivativePtr<CompositeVector>(kr_key_, tag_next_, temp_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> kr_uw =
        S_->GetPtr<CompositeVector>(kr_uw_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> flux =
        S_->GetPtr<CompositeVector>(water_flux_key_, tag_next_);

      S_->GetEvaluator(potential_key_, tag_next_).Update(*S_, name_);
      Teuchos::RCP<const CompositeVector> pres_elev =
        S_->GetPtr<CompositeVector>(potential_key_, tag_next_);

      // form the operator
      ddivq_dT_->SetScalarCoefficient(kr_uw, dkrdT);
      ddivq_dT_->UpdateMatrices(flux.ptr(), pres_elev.ptr());
      ddivq_dT_->UpdateMatricesNewtonCorrection(flux.ptr(), pres_elev.ptr());
      ddivq_dT_->ApplyBCs(false, true, false);
    }

    // -- dE/dp diagonal term
    S_->GetEvaluator(e_key_, tag_next_).UpdateDerivative(*S_, name_, pres_key_, tag_next_);
    Teuchos::RCP<const CompositeVector> dE_dp =
      S_->GetDerivativePtr<CompositeVector>(e_key_, tag_next_, pres_key_, tag_next_);

    // -- scale to dE/dh
    S_->GetEvaluator(pd_bar_key_, tag_next_).UpdateDerivative(*S_, name_, pres_key_, tag_next_);
    auto dh_dp =
      S_->GetDerivativePtr<CompositeVector>(pd_bar_key_, tag_next_, pres_key_, tag_next_);

    // -- add it in
    CompositeVector dE_dh(dE_dp->Map());
    dE_dh.ReciprocalMultiply(1. / h, *dh_dp, *dE_dp, 0.);
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
int
MPCSurface::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Precon application:" << std::endl;

  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Residuals:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  r_p");
    vnames.push_back("  r_T");
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(u->SubVector(0)->Data().ptr());
    vecs.push_back(u->SubVector(1)->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  int ierr = 0;
  if (precon_type_ == PRECON_NONE) {
    *Pu = *u;
    ierr = 1;
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    ierr = StrongMPC::ApplyPreconditioner(u, Pu);
  } else if (precon_type_ == PRECON_PICARD || (precon_type_ == PRECON_EWC)) {
    ierr = preconditioner_->ApplyInverse(*u, *Pu);

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "PC * residuals:" << std::endl;
      std::vector<std::string> vnames;
      vnames.push_back("  PC*r_h");
      vnames.push_back("  PC*r_T");
      std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
      vecs.push_back(Pu->SubVector(0)->Data().ptr());
      vecs.push_back(Pu->SubVector(1)->Data().ptr());
      db_->WriteVectors(vnames, vecs, true);
    }

    // tack on the variable change from h to p
    const Epetra_MultiVector& dh_dp =
      *S_->GetDerivativePtr<CompositeVector>(pd_bar_key_, tag_next_, pres_key_, tag_next_)
         ->ViewComponent("cell", false);
    Epetra_MultiVector& Pu_c = *Pu->SubVector(0)->Data()->ViewComponent("cell", false);

    for (unsigned int c = 0; c != Pu_c.MyLength(); ++c) { Pu_c[0][c] /= dh_dp[0][c]; }
  }

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "PC * residuals:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  PC*r_p");
    vnames.push_back("  PC*r_T");
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vecs.push_back(Pu->SubVector(0)->Data().ptr());
    vecs.push_back(Pu->SubVector(1)->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  return (ierr > 0) ? 0 : 1;
}


} // namespace Amanzi
