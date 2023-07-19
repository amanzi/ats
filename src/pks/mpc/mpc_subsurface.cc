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

#include "permafrost_model.hh"
#include "liquid_ice_model.hh"
#include "richards.hh"
#include "mpc_delegate_ewc_subsurface.hh"
#include "mpc_subsurface.hh"

#define DEBUG_FLAG 1

namespace Amanzi {

MPCSubsurface::MPCSubsurface(Teuchos::ParameterList& pk_tree_list,
                             const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree_list, global_list, S, soln),
    StrongMPC<PK_PhysicalBDF_Default>(pk_tree_list, global_list, S, soln),
    update_pcs_(0)
{
  // set up keys
  dump_ = plist_->get<bool>("dump preconditioner", false);
  domain_name_ = plist_->get<std::string>("domain name", "domain");

  temp_key_ = Keys::readKey(*plist_, domain_name_, "temperature", "temperature");
  pres_key_ = Keys::readKey(*plist_, domain_name_, "pressure", "pressure");
  e_key_ = Keys::readKey(*plist_, domain_name_, "energy", "energy");
  wc_key_ = Keys::readKey(*plist_, domain_name_, "water content", "water_content");
  tc_key_ = Keys::readKey(*plist_, domain_name_, "thermal conductivity", "thermal_conductivity");
  uw_tc_key_ = Keys::readKey(
    *plist_, domain_name_, "upwinded thermal conductivity", "upwind_thermal_conductivity");
  kr_key_ = Keys::readKey(*plist_, domain_name_, "conductivity", "relative_permeability");
  uw_kr_key_ =
    Keys::readKey(*plist_, domain_name_, "upwinded conductivity", "upwind_relative_permeability");
  enth_key_ = Keys::readKey(*plist_, domain_name_, "enthalpy", "enthalpy");
  hkr_key_ = Keys::readKey(
    *plist_, domain_name_, "enthalpy times conductivity", "enthalpy_times_relative_permeability");
  uw_hkr_key_ = Keys::readKey(*plist_,
                              domain_name_,
                              "upwind_enthalpy times conductivity",
                              "upwind_enthalpy_times_relative_permeability");
  energy_flux_key_ =
    Keys::readKey(*plist_, domain_name_, "diffusive energy flux", "diffusive_energy_flux");
  water_flux_key_ = Keys::readKey(*plist_, domain_name_, "water flux", "water_flux");
  water_flux_dir_key_ =
    Keys::readKey(*plist_, domain_name_, "water flux direction", "water_flux_direction");
  rho_key_ = Keys::readKey(*plist_, domain_name_, "mass density liquid", "mass_density_liquid");
}

// -- Initialize owned (dependent) variables.
void
MPCSubsurface::Setup()
{
  auto pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
  global_list_->sublist("PKs").sublist(pk_order[0]).set("scale preconditioner to pressure", false);

  // supress energy's vision of advective terms as we can do better
  if (!plist_->get<bool>("supress Jacobian terms: d div hq / dp,T", false)) {
    if (pks_list_->sublist(pk_order[1]).isParameter("supress advective terms in preconditioner") &&
        !pks_list_->sublist(pk_order[1]).get("supress advective terms in preconditioner", false)) {
      Errors::Message msg("MPC Incorrect input: options \"supress Jacobian terms: d div hq / "
                          "dp,T\" and subsurface energy PK option \"supress advective terms in "
                          "preconditioner\" should not both be false, as these include some of the "
                          "same Jacobian information.\n Recommended: Enable/suppress the latter.");

      Exceptions::amanzi_throw(msg);
    }
  }

  // set up the sub-pks
  StrongMPC<PK_PhysicalBDF_Default>::Setup();
  mesh_ = S_->GetMesh(domain_name_);

  // set up debugger
  db_ = sub_pks_[0]->debugger();

  S_->Require<CompositeVector, CompositeVectorSpace>(rho_key_, tag_next_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(rho_key_, tag_next_);

  // see amanzi/ats#167
  // if (S_->GetEvaluator(e_key_, tag_next_).IsDifferentiableWRT(*S_, pres_key_, tag_next_)) {
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    e_key_, tag_next_, pres_key_, tag_next_, e_key_);
  // }
  // if (S_->GetEvaluator(wc_key_, tag_next_).IsDifferentiableWRT(*S_, temp_key_, tag_next_)) {
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    wc_key_, tag_next_, temp_key_, tag_next_, wc_key_);
  // }

  // Get the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Operator> pcA = sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Operator> pcB = sub_pks_[1]->preconditioner();

  if (pks_list_->sublist(pk_order[0]).isSublist("diffusion") &&
      pks_list_->sublist(pk_order[0])
          .sublist("diffusion")
          .get<std::string>("discretization primary") == "fv: default" &&
      pks_list_->sublist(pk_order[1]).isSublist("diffusion") &&
      pks_list_->sublist(pk_order[1])
          .sublist("diffusion")
          .get<std::string>("discretization primary") == "fv: default") {
    is_fv_ = true;
  } else {
    is_fv_ = false;
  }

  // Create the combined operator
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcA->DomainMap()))));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcB->DomainMap()))));

  preconditioner_ = Teuchos::rcp(new Operators::TreeOperator(tvs, plist_->sublist("operator")));
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
    std::vector<AmanziMesh::Entity_kind> locations2(2);
    std::vector<std::string> names2(2);
    std::vector<int> num_dofs2(2, 1);
    locations2[0] = AmanziMesh::Entity_kind::CELL;
    names2[0] = "cell";

    // Create the block for derivatives of mass conservation with respect to temperature
    // -- derivatives of kr with respect to temperature
    if (precon_type_ != PRECON_NO_FLOW_COUPLING &&
        !plist_->get<bool>("supress Jacobian terms: d div q / dT", false)) {
      // require the derivative drel_perm/dT
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        kr_key_, tag_next_, temp_key_, tag_next_);

      if (!is_fv_) {
        duw_krdT_key_ = Keys::getDerivKey(uw_kr_key_, temp_key_);
        S_->Require<CompositeVector, CompositeVectorSpace>(duw_krdT_key_, tag_next_, name_)
          .SetMesh(mesh_)
          ->SetGhosted()
          ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
        S_->GetRecordW(duw_krdT_key_, tag_next_, name_).set_io_vis(false);

        upwinding_dkrdT_ = Teuchos::rcp(
          new Operators::UpwindTotalFlux(name_, tag_next_, water_flux_dir_key_, 1.e-8));
      }

      // set up the operator
      Teuchos::ParameterList divq_plist(
        pks_list_->sublist(pk_order[0]).sublist("diffusion preconditioner"));

      if (is_fv_)
        divq_plist.set("Newton correction", "true Jacobian");
      else
        divq_plist.set("Newton correction", "approximate Jacobian");

      divq_plist.set("exclude primary terms", true);

      Operators::PDE_DiffusionFactory opfactory;

      ddivq_dT_ = opfactory.CreateWithGravity(divq_plist, mesh_);
      dWC_dT_block_ = ddivq_dT_->global_operator();
    }

    // -- derivatives of water content with respect to temperature
    // see amanzi/ats#167
    // if (S_->GetEvaluator(wc_key_, tag_next_).IsDifferentiableWRT(*S_, temp_key_, tag_next_)) {
    if (dWC_dT_block_ == Teuchos::null) {
      dWC_dT_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
      dWC_dT_block_ = dWC_dT_->global_operator();
    } else {
      dWC_dT_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, dWC_dT_block_));
    }
    // }

    // Create the block for derivatives of energy conservation with respect to pressure
    // -- derivatives of thermal conductivity with respect to pressure
    if (precon_type_ != PRECON_NO_FLOW_COUPLING &&
        !plist_->get<bool>("supress Jacobian terms: d div K grad T / dp", false)) {
      // require the derivative dkappa/dp
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        tc_key_, tag_next_, pres_key_, tag_next_);

      // need to upwind dkappa/dp
      if (!is_fv_) {
        duw_tcdp_key_ = Keys::getDerivKey(uw_tc_key_, pres_key_);
        S_->Require<CompositeVector, CompositeVectorSpace>(duw_tcdp_key_, tag_next_, name_)
          .SetMesh(mesh_)
          ->SetGhosted()
          ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
        S_->GetRecordW(duw_tcdp_key_, tag_next_, name_).set_io_vis(false);

        upwinding_dkappa_dp_ =
          Teuchos::rcp(new Operators::UpwindTotalFlux(name_, tag_next_, energy_flux_key_, 1.e-8));
      }

      // set up the operator
      Teuchos::ParameterList ddivKgT_dp_plist(
        pks_list_->sublist(pk_order[1]).sublist("diffusion preconditioner"));
      if (is_fv_)
        ddivKgT_dp_plist.set("Newton correction", "true Jacobian");
      else
        ddivKgT_dp_plist.set("Newton correction", "approximate Jacobian");

      ddivKgT_dp_plist.set("exclude primary terms", true);

      Operators::PDE_DiffusionFactory opfactory;

      if (dE_dp_block_ == Teuchos::null) {
        ddivKgT_dp_ = opfactory.Create(ddivKgT_dp_plist, mesh_);
        dE_dp_block_ = ddivKgT_dp_->global_operator();
      } else {
        ddivKgT_dp_ = opfactory.Create(ddivKgT_dp_plist, dE_dp_block_);
      }
      ddivKgT_dp_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[0]->BCs());
    }


    // -- derivatives of advection term
    if (precon_type_ != PRECON_NO_FLOW_COUPLING &&
        !plist_->get<bool>("supress Jacobian terms: d div hq / dp,T", false)) {
      // derivative with respect to pressure
      Teuchos::ParameterList divhq_dp_plist(
        pks_list_->sublist(pk_order[0]).sublist("diffusion preconditioner"));

      if (is_fv_)
        divhq_dp_plist.set("Newton correction", "true Jacobian");
      else
        divhq_dp_plist.set("Newton correction", "approximate Jacobian");

      Operators::PDE_DiffusionFactory opfactory;
      if (dE_dp_block_ == Teuchos::null) {
        ddivhq_dp_ = opfactory.CreateWithGravity(divhq_dp_plist, mesh_);
        dE_dp_block_ = ddivhq_dp_->global_operator();
      } else {
        ddivhq_dp_ = opfactory.CreateWithGravity(divhq_dp_plist, dE_dp_block_);
      }

      // derivative with respect to temperature
      Teuchos::ParameterList divhq_dT_plist(
        pks_list_->sublist(pk_order[0]).sublist("diffusion preconditioner"));
      divhq_dT_plist.set("exclude primary terms", true);

      if (is_fv_)
        divhq_dT_plist.set("Newton correction", "true Jacobian");
      else
        divhq_dT_plist.set("Newton correction", "approximate Jacobian");

      ddivhq_dT_ = opfactory.CreateWithGravity(divhq_dT_plist, pcB);

      // need a field, evaluator, and upwinding for h * kr * rho/mu
      // -- first the evaluator
      auto& hkr_eval_list = S_->GetEvaluatorList(hkr_key_);
      hkr_eval_list.set("evaluator name", hkr_key_);
      Teuchos::Array<std::string> deps(2);
      deps[0] = enth_key_;
      deps[1] = kr_key_;
      hkr_eval_list.set("dependencies", deps);
      hkr_eval_list.set("evaluator type", "multiplicative evaluator");

      // -- now the field
      names2[1] = "boundary_face";
      locations2[1] = AmanziMesh::Entity_kind::BOUNDARY_FACE;
      S_->Require<CompositeVector, CompositeVectorSpace>(hkr_key_, tag_next_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->AddComponents(names2, locations2, num_dofs2);
      S_->RequireEvaluator(hkr_key_, tag_next_);

      S_->Require<CompositeVector, CompositeVectorSpace>(uw_hkr_key_, tag_next_, name_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
      S_->GetRecordW(uw_hkr_key_, tag_next_, name_).set_io_vis(false);

      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        hkr_key_, tag_next_, temp_key_, tag_next_);
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        hkr_key_, tag_next_, pres_key_, tag_next_);

      std::string method_name =
        pks_list_->sublist(pk_order[0])
          .get<std::string>("relative permeability method", "upwind with gravity");
      if (method_name != "upwind with Darcy flux") {
        Errors::Message msg;
        msg << "Subsurface coupler with advective Jacobian terms only supports a Richards upwind "
               "scheme of "
            << "\"upwind with Darcy flux\", but the method \"" << method_name
            << "\" was requested.";
        Exceptions::amanzi_throw(msg);
      }
      upwinding_hkr_ =
        Teuchos::rcp(new Operators::UpwindTotalFlux(name_, tag_next_, water_flux_dir_key_, 1.e-8));

      if (!is_fv_) {
        // -- and the upwinded field
        locations2[1] = AmanziMesh::Entity_kind::FACE;
        names2[1] = "face";

        S_->Require<CompositeVector, CompositeVectorSpace>(
            Keys::getDerivKey(uw_hkr_key_, pres_key_), tag_next_, name_)
          .SetMesh(mesh_)
          ->SetGhosted()
          ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
        S_->GetRecordW(Keys::getDerivKey(uw_hkr_key_, pres_key_), tag_next_, name_)
          .set_io_vis(false);

        S_->Require<CompositeVector, CompositeVectorSpace>(
            Keys::getDerivKey(uw_hkr_key_, temp_key_), tag_next_, name_)
          .SetMesh(mesh_)
          ->SetGhosted()
          ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
        S_->GetRecordW(Keys::getDerivKey(uw_hkr_key_, temp_key_), tag_next_, name_)
          .set_io_vis(false);

        // -- and the upwinding
        upwinding_dhkr_dp_ = Teuchos::rcp(
          new Operators::UpwindTotalFlux(name_, tag_next_, water_flux_dir_key_, 1.e-8));
        upwinding_dhkr_dT_ = Teuchos::rcp(
          new Operators::UpwindTotalFlux(name_, tag_next_, water_flux_dir_key_, 1.e-8));
      }
    }

    // -- derivatives of energy with respect to pressure
    // see amanzi/ats#167
    // if (S_->GetEvaluator(wc_key_, tag_next_).IsDifferentiableWRT(*S_, temp_key_, tag_next_)) {
    if (dE_dp_block_ == Teuchos::null) {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
      dE_dp_block_ = dE_dp_->global_operator();
    } else {
      dE_dp_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, dE_dp_block_));
    }
    // }

    AMANZI_ASSERT(dWC_dT_block_ != Teuchos::null);
    AMANZI_ASSERT(dE_dp_block_ != Teuchos::null);
    preconditioner_->set_operator_block(0, 1, dWC_dT_block_);
    preconditioner_->set_operator_block(1, 0, dE_dp_block_);

    // set up sparsity structure
    preconditioner_->set_inverse_parameters(plist_->sublist("inverse"));
  }

  // create the EWC delegate
  if (plist_->isSublist("ewc delegate")) {
    Teuchos::RCP<Teuchos::ParameterList> sub_ewc_list = Teuchos::sublist(plist_, "ewc delegate");
    sub_ewc_list->set("PK name", name_);
    sub_ewc_list->set("domain name", domain_name_);
    ewc_ = Teuchos::rcp(new MPCDelegateEWCSubsurface(*sub_ewc_list, S_));
    ewc_->set_tags(tag_current_, tag_next_);

    Teuchos::RCP<EWCModelBase> model;
    if (S_->HasRecordSet("internal_energy_gas")) {
      model = Teuchos::rcp(new PermafrostModel());
    } else {
      model = Teuchos::rcp(new LiquidIceModel());
    }
    ewc_->set_model(model);
    ewc_->setup();
  }
}


void
MPCSubsurface::Initialize()
{
  StrongMPC<PK_PhysicalBDF_Default>::Initialize();
  if (ewc_ != Teuchos::null) ewc_->initialize();

  // initialize offdiagonal operators
  if (precon_type_ != PRECON_NONE && precon_type_ != PRECON_BLOCK_DIAGONAL) {
    if (S_->HasDerivative(wc_key_, tag_next_, temp_key_, tag_next_))
      S_->GetDerivativeW<CompositeVector>(wc_key_, tag_next_, temp_key_, tag_next_, wc_key_)
        .PutScalar(0.0);

    if (S_->HasDerivative(e_key_, tag_next_, pres_key_, tag_next_))
      S_->GetDerivativeW<CompositeVector>(e_key_, tag_next_, pres_key_, tag_next_, e_key_)
        .PutScalar(0.0);
  }

  Teuchos::RCP<Flow::Richards> richards_pk;
  if (ddivq_dT_ != Teuchos::null) {
    if (richards_pk == Teuchos::null) {
      richards_pk = Teuchos::rcp_dynamic_cast<Flow::Richards>(sub_pks_[0]);
      AMANZI_ASSERT(richards_pk != Teuchos::null);
    }

    if (!is_fv_) {
      S_->GetW<CompositeVector>(duw_krdT_key_, tag_next_, name_).PutScalar(0.0);
      S_->GetRecordW(duw_krdT_key_, tag_next_, name_).set_initialized();
    }

    const AmanziGeometry::Point& g = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    ddivq_dT_->SetGravity(g);
    ddivq_dT_->SetBCs(sub_pks_[0]->BCs(), sub_pks_[1]->BCs());
    ddivq_dT_->SetTensorCoefficient(richards_pk->K_);
  }

  if (ddivKgT_dp_ != Teuchos::null) {
    if (!is_fv_) {
      S_->GetW<CompositeVector>(duw_tcdp_key_, tag_next_, name_).PutScalar(0.0);
      S_->GetRecordW(duw_tcdp_key_, tag_next_, name_).set_initialized();
    }

    ddivKgT_dp_->SetTensorCoefficient(Teuchos::null);
  }

  if (ddivhq_dp_ != Teuchos::null) {
    if (richards_pk == Teuchos::null) {
      richards_pk = Teuchos::rcp_dynamic_cast<Flow::Richards>(sub_pks_[0]);
      AMANZI_ASSERT(richards_pk != Teuchos::null);
    }

    S_->GetW<CompositeVector>(uw_hkr_key_, tag_next_, name_).PutScalar(1.0);
    S_->GetRecordW(uw_hkr_key_, tag_next_, name_).set_initialized();

    if (!is_fv_) {
      S_->GetW<CompositeVector>(Keys::getDerivKey(uw_hkr_key_, pres_key_), tag_next_, name_)
        .PutScalar(0.0);
      S_->GetRecordW(Keys::getDerivKey(uw_hkr_key_, pres_key_), tag_next_, name_).set_initialized();
      S_->GetW<CompositeVector>(Keys::getDerivKey(uw_hkr_key_, temp_key_), tag_next_, name_)
        .PutScalar(0.0);
      S_->GetRecordW(Keys::getDerivKey(uw_hkr_key_, temp_key_), tag_next_, name_).set_initialized();
    }

    const AmanziGeometry::Point& g = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    ddivhq_dp_->SetGravity(g);
    ddivhq_dp_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[0]->BCs());
    ddivhq_dp_->SetTensorCoefficient(richards_pk->K_);

    ddivhq_dT_->SetGravity(g);
    ddivhq_dT_->SetBCs(sub_pks_[1]->BCs(), sub_pks_[1]->BCs());
    ddivhq_dT_->SetTensorCoefficient(richards_pk->K_);
  }
}


void
MPCSubsurface::set_tags(const Tag& tag_current, const Tag& tag_next)
{
  StrongMPC<PK_PhysicalBDF_Default>::set_tags(tag_current, tag_next);
  if (ewc_ != Teuchos::null) ewc_->set_tags(tag_current, tag_next);
}


void
MPCSubsurface::CommitStep(double t_old, double t_new, const Tag& tag)
{
  if (ewc_ != Teuchos::null) { ewc_->commit_state(); }

  update_pcs_ = 0;
  StrongMPC<PK_PhysicalBDF_Default>::CommitStep(t_old, t_new, tag);
}


// update the predictor to be physically consistent
bool
MPCSubsurface::ModifyPredictor(double h,
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
MPCSubsurface::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  if (precon_type_ == PRECON_NONE) {
    // nothing to do
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    StrongMPC::UpdatePreconditioner(t, up, h);
  } else if (precon_type_ == PRECON_PICARD || precon_type_ == PRECON_EWC) {
    preconditioner_->InitOffdiagonals(); // zero out offdiagonal blocks and mark for re-computation
    StrongMPC::UpdatePreconditioner(t, up, h);

    // dWC / dT block
    // -- dkr/dT
    // -- update and upwind d kr / dT
    if (ddivq_dT_ != Teuchos::null) {
      S_->GetEvaluator(kr_key_, tag_next_).UpdateDerivative(*S_, name_, temp_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> dkrdT =
        S_->GetDerivativePtr<CompositeVector>(kr_key_, tag_next_, temp_key_, tag_next_);
      if (!is_fv_) {
        Teuchos::RCP<CompositeVector> duw_kr =
          S_->GetPtrW<CompositeVector>(duw_krdT_key_, tag_next_, name_);
        duw_kr->PutScalar(0.0);
        upwinding_dkrdT_->Update(*dkrdT, *duw_kr, *S_);
        dkrdT = S_->GetPtr<CompositeVector>(duw_krdT_key_, tag_next_);
      }

      // form the operator
      Teuchos::RCP<const CompositeVector> kr_uw =
        S_->GetPtr<CompositeVector>(uw_kr_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> flux =
        S_->GetPtr<CompositeVector>(water_flux_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(rho_key_, tag_next_);

      ddivq_dT_->SetDensity(rho);
      ddivq_dT_->SetScalarCoefficient(kr_uw, dkrdT);
      ddivq_dT_->UpdateMatrices(flux.ptr(), up->SubVector(0)->Data().ptr());
      ddivq_dT_->UpdateMatricesNewtonCorrection(flux.ptr(), up->SubVector(0)->Data().ptr());

      ddivq_dT_->ApplyBCs(false, true, false);
    }

    // -- dWC/dT diagonal term
    Teuchos::RCP<const CompositeVector> dWC_dT = Teuchos::null;
    if (S_->GetEvaluator(wc_key_, tag_next_).IsDifferentiableWRT(*S_, temp_key_, tag_next_)) {
      S_->GetEvaluator(wc_key_, tag_next_).UpdateDerivative(*S_, name_, temp_key_, tag_next_);
      dWC_dT = S_->GetDerivativePtr<CompositeVector>(wc_key_, tag_next_, temp_key_, tag_next_);
      dWC_dT_->AddAccumulationTerm(*dWC_dT, h, "cell", false);
    }

    // dE / dp block
    // -- d Kappa / dp
    // Update and upwind thermal conductivity
    if (ddivKgT_dp_ != Teuchos::null) {
      S_->GetEvaluator(tc_key_, tag_next_).UpdateDerivative(*S_, name_, pres_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> dkappa_dp =
        S_->GetDerivativePtr<CompositeVector>(tc_key_, tag_next_, pres_key_, tag_next_);
      if (!is_fv_) {
        Teuchos::RCP<CompositeVector> duw_kappa_dp =
          S_->GetPtrW<CompositeVector>(duw_tcdp_key_, tag_next_, name_);
        duw_kappa_dp->PutScalar(0.0);
        upwinding_dkappa_dp_->Update(*dkappa_dp, *duw_kappa_dp, *S_);
        dkappa_dp = S_->GetPtr<CompositeVector>(duw_tcdp_key_, tag_next_);
      }

      // form the operator
      Teuchos::RCP<const CompositeVector> uw_Kappa =
        S_->GetPtr<CompositeVector>(uw_tc_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> flux =
        S_->GetPtr<CompositeVector>(energy_flux_key_, tag_next_);
      ddivKgT_dp_->SetScalarCoefficient(uw_Kappa, dkappa_dp);
      ddivKgT_dp_->UpdateMatrices(flux.ptr(), up->SubVector(1)->Data().ptr());
      ddivKgT_dp_->UpdateMatricesNewtonCorrection(flux.ptr(), up->SubVector(1)->Data().ptr());

      ddivKgT_dp_->ApplyBCs(false, true, false);
    }

    // -- d adv / dp   This one is a bit more complicated...
    // Update and upwind enthalpy * kr * rho/mu
    if (ddivhq_dp_ != Teuchos::null) {
      // -- update values
      S_->GetEvaluator(hkr_key_, tag_next_).Update(*S_, name_);
      S_->GetEvaluator(hkr_key_, tag_next_).UpdateDerivative(*S_, name_, pres_key_, tag_next_);
      S_->GetEvaluator(hkr_key_, tag_next_).UpdateDerivative(*S_, name_, temp_key_, tag_next_);

      Teuchos::RCP<const CompositeVector> denth_kr_dp_uw;
      Teuchos::RCP<const CompositeVector> denth_kr_dT_uw;

      Teuchos::RCP<const CompositeVector> enth_kr =
        S_->GetPtr<CompositeVector>(hkr_key_, tag_next_);
      Teuchos::RCP<CompositeVector> enth_kr_uw =
        S_->GetPtrW<CompositeVector>(uw_hkr_key_, tag_next_, name_);
      enth_kr_uw->PutScalar(0.0);

      enth_kr_uw->ViewComponent("face", false)
        ->Export(
          *enth_kr->ViewComponent("boundary_face", false), mesh_->getBoundaryFaceImporter(), Insert);
      upwinding_hkr_->Update(*enth_kr, *enth_kr_uw, *S_);

      // -- stick zeros in the boundary faces
      Epetra_MultiVector enth_kr_bf(*enth_kr->ViewComponent("boundary_face", false));
      enth_kr_bf.PutScalar(0.0);
      enth_kr_uw->ViewComponent("face", false)
        ->Export(enth_kr_bf, mesh_->getBoundaryFaceImporter(), Insert);

      if (is_fv_) {
        denth_kr_dp_uw =
          S_->GetDerivativePtr<CompositeVector>(hkr_key_, tag_next_, pres_key_, tag_next_);
        denth_kr_dT_uw =
          S_->GetDerivativePtr<CompositeVector>(hkr_key_, tag_next_, temp_key_, tag_next_);
      } else {
        Teuchos::RCP<const CompositeVector> denth_kr_dp =
          S_->GetDerivativePtr<CompositeVector>(hkr_key_, tag_next_, pres_key_, tag_next_);
        Teuchos::RCP<const CompositeVector> denth_kr_dT =
          S_->GetDerivativePtr<CompositeVector>(hkr_key_, tag_next_, temp_key_, tag_next_);

        // -- zero target data (may be unnecessary?)
        Teuchos::RCP<CompositeVector> denth_kr_dp_uw_nc =
          S_->GetPtrW<CompositeVector>(Keys::getDerivKey(uw_hkr_key_, pres_key_), tag_next_, name_);
        denth_kr_dp_uw_nc->PutScalar(0.);
        Teuchos::RCP<CompositeVector> denth_kr_dT_uw_nc =
          S_->GetPtrW<CompositeVector>(Keys::getDerivKey(uw_hkr_key_, temp_key_), tag_next_, name_);
        denth_kr_dT_uw_nc->PutScalar(0.);

        // -- copy boundary faces into upwinded vector
        denth_kr_dp_uw_nc->ViewComponent("face", false)
          ->Export(*denth_kr_dp->ViewComponent("boundary_face", false),
                   mesh_->getBoundaryFaceImporter(),
                   Insert);
        denth_kr_dT_uw_nc->ViewComponent("face", false)
          ->Export(*denth_kr_dT->ViewComponent("boundary_face", false),
                   mesh_->getBoundaryFaceImporter(),
                   Insert);

        // -- upwind
        upwinding_dhkr_dp_->Update(*denth_kr_dp, *denth_kr_dp_uw_nc, *S_);
        upwinding_dhkr_dT_->Update(*denth_kr_dT, *denth_kr_dT_uw_nc, *S_);

        // -- stick zeros in the boundary faces
        Epetra_MultiVector enth_kr_bf(*enth_kr->ViewComponent("boundary_face", false));
        enth_kr_bf.PutScalar(0.0);
        denth_kr_dp_uw_nc->ViewComponent("face", false)
          ->Export(enth_kr_bf, mesh_->getBoundaryFaceImporter(), Insert);
        denth_kr_dT_uw_nc->ViewComponent("face", false)
          ->Export(enth_kr_bf, mesh_->getBoundaryFaceImporter(), Insert);

        denth_kr_dp_uw =
          S_->GetPtr<CompositeVector>(Keys::getDerivKey(uw_hkr_key_, pres_key_), tag_next_);
        denth_kr_dT_uw =
          S_->GetPtr<CompositeVector>(Keys::getDerivKey(uw_hkr_key_, temp_key_), tag_next_);
      }

      Teuchos::RCP<const CompositeVector> flux =
        S_->GetPtr<CompositeVector>(water_flux_key_, tag_next_);
      Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(rho_key_, tag_next_);

      // form the operator: pressure component
      ddivhq_dp_->SetDensity(rho);
      ddivhq_dp_->SetScalarCoefficient(enth_kr_uw, denth_kr_dp_uw);
      // -- update the local matrices, div h * kr grad
      ddivhq_dp_->UpdateMatrices(Teuchos::null, Teuchos::null);
      // -- determine the advective fluxes, q_a = h * kr grad p
      CompositeVector adv_flux(*flux, INIT_MODE_ZERO);
      Teuchos::Ptr<CompositeVector> adv_flux_ptr(&adv_flux);
      ddivhq_dp_->UpdateFlux(up->SubVector(0)->Data().ptr(), adv_flux_ptr);
      // -- add in components div (d h*kr / dp) grad q_a / (h*kr)
      ddivhq_dp_->UpdateMatricesNewtonCorrection(adv_flux_ptr, up->SubVector(0)->Data().ptr());
      ddivhq_dp_->ApplyBCs(false, true, false);

      // form the operator: temperature component
      ddivhq_dT_->SetDensity(rho);
      ddivhq_dT_->SetScalarCoefficient(enth_kr_uw, denth_kr_dT_uw);
      // -- add in components div (d h*kr / dp) grad q_a / (h*kr)
      ddivhq_dT_->UpdateMatrices(adv_flux_ptr, up->SubVector(0)->Data().ptr());
      ddivhq_dT_->UpdateMatricesNewtonCorrection(adv_flux_ptr, up->SubVector(0)->Data().ptr());
      ddivhq_dT_->ApplyBCs(false, true, false);
    }

    // -- dE/dp diagonal term
    Teuchos::RCP<const CompositeVector> dE_dp = Teuchos::null;
    if (S_->GetEvaluator(e_key_, tag_next_).IsDifferentiableWRT(*S_, pres_key_, tag_next_)) {
      S_->GetEvaluator(e_key_, tag_next_).UpdateDerivative(*S_, name_, pres_key_, tag_next_);
      dE_dp = S_->GetDerivativePtr<CompositeVector>(e_key_, tag_next_, pres_key_, tag_next_);
      dE_dp_->AddAccumulationTerm(*dE_dp, h, "cell", false);
    }

    // write for debugging
    std::vector<std::string> vnames;
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    if (dWC_dT != Teuchos::null) {
      vnames.push_back("  dwc_dT");
      vecs.push_back(dWC_dT.ptr());
    }
    if (dE_dp != Teuchos::null) {
      vnames.push_back("  de_dp");
      vecs.push_back(dE_dp.ptr());
    }
    if (vecs.size() > 0) db_->WriteVectors(vnames, vecs, false);
  }

  if (precon_type_ == PRECON_EWC) { ewc_->UpdatePreconditioner(t, up, h); }
  update_pcs_++;
}


// -----------------------------------------------------------------------------
// Wrapper to call the requested preconditioner.
// -----------------------------------------------------------------------------
int
MPCSubsurface::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
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
  } else if (precon_type_ == PRECON_PICARD) {
    ierr = preconditioner_->ApplyInverse(*u, *Pu);
  } else if (precon_type_ == PRECON_EWC) {
    ierr = preconditioner_->ApplyInverse(*u, *Pu);
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
