/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "pk_helpers.hh"

#include "surface_balance_base.hh"

namespace Amanzi {
namespace SurfaceBalance {


SurfaceBalanceBase::SurfaceBalanceBase(Teuchos::ParameterList& pk_tree,
                                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                       const Teuchos::RCP<State>& S,
                                       const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, global_list, S, solution), PK_PhysicalBDF_Default(pk_tree, global_list, S, solution)
{
  // source terms
  eps_ = plist_->get<double>("source term finite difference epsilon", 1.e-8);
  is_source_ = plist_->get<bool>("source term", true);
  if (is_source_ && source_key_.empty()) {
    source_key_ = Keys::readKey(*plist_, domain_, "source", "source_sink");
  }
  source_finite_difference_ = plist_->get<bool>("source term finite difference", false);
  is_source_differentiable_ = plist_->get<bool>("source term is differentiable", true);

  modify_predictor_positivity_preserving_ =
    plist_->get<bool>("modify predictor positivity preserving", false);

  theta_ = plist_->get<double>("time discretization theta", 1.0);
  if (theta_ > 1 || theta_ < 0) {
    Errors::Message message(
      "SurfaceBalanceBase: \"time discretization theta\" value must be between 0 and 1.");
    Exceptions::amanzi_throw(message);
  }

  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .01 * 55000.); // h * nl
}


// main methods
// -- Setup data.
void
SurfaceBalanceBase::Setup()
{
  PK_PhysicalBDF_Default::Setup();

  // requirements: primary variable
  //  NOTE: no need to require evaluator here, either at the old or new
  //  times, as this was done in pk_physical.  All we have to do is set the
  //  structure.
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirements: source terms from above
  if (is_source_) {
    if (theta_ > 0) {
      requireAtNext(source_key_, tag_next_, *S_)
        .SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

      if (is_source_differentiable_ && !source_finite_difference_
          // this cannot work in general yet, see amanzi/ats#167
          //&& S_->GetEvaluator(source_key_, tag_next_).IsDifferentiableWRT(*S_, key_, tag_next_)
      ) {
        S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
            source_key_, tag_next_, key_, tag_next_)
          .SetMesh(mesh_)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
        // NOTE, remove SetMesh/AddComponent lines after fixing amanzi/ats#167.
        // The mesh should get set by the evaluator, but when
        // the evaluator isn't actually differentiable, it
        // doesn't get done.
      }
    }
    if (theta_ < 1) {
      requireAtCurrent(source_key_, tag_current_, *S_, name_)
        .SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

  // requirements: conserved quantity at current and new times
  conserved_quantity_ = conserved_key_ != key_;
  if (conserved_quantity_) {
    requireAtNext(conserved_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
      conserved_key_, tag_next_, key_, tag_next_);

    //    and at the current time, where it is a copy evaluator
    requireAtCurrent(conserved_key_, tag_current_, *S_, name_);
  }

  // operator for inverse
  Teuchos::ParameterList& acc_plist = plist_->sublist("accumulation preconditioner");
  acc_plist.set("entity kind", "cell");
  acc_plist.set("inverse", plist_->sublist("inverse"));
  // old style... deprecate me!
  acc_plist.sublist("inverse").setParameters(plist_->sublist("preconditioner"));
  acc_plist.sublist("inverse").setParameters(plist_->sublist("linear solver"));

  preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_plist, mesh_));
  preconditioner_ = preconditioner_acc_->global_operator();
}


void
SurfaceBalanceBase::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  // saves primary variable
  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, tag_next);

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;

  if (theta_ < 1.0) assign(source_key_, tag_current, tag_next, *S_);
}


// computes the non-linear functional g = g(t,u,udot)
void
SurfaceBalanceBase::FunctionalResidual(double t_old,
                                       double t_new,
                                       Teuchos::RCP<TreeVector> u_old,
                                       Teuchos::RCP<TreeVector> u_new,
                                       Teuchos::RCP<TreeVector> g)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  double dt = t_new - t_old;

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, tag_next_);

  bool debug = false;
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) debug = true;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old << " t1 = " << t_new << " h = " << dt
               << std::endl;
    std::vector<std::string> vnames;
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vnames.push_back("u_old");
    vnames.push_back("u_new");
    vecs.push_back(S_->GetPtr<CompositeVector>(key_, tag_current_).ptr());
    vecs.push_back(u_new->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  if (conserved_quantity_) {
    S_->GetEvaluator(conserved_key_, tag_next_).Update(*S_, name_);
    // S_->GetEvaluator(conserved_key_, tag_current_).Update(*S_, name_);
    auto& conserved1 = S_->Get<CompositeVector>(conserved_key_, tag_next_);
    auto& conserved0 = S_->Get<CompositeVector>(conserved_key_, tag_current_);
    g->Data()->Update(1.0 / dt, conserved1, -1.0 / dt, conserved0, 0.0);

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      std::vector<std::string> vnames;
      std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
      vnames.push_back("C_old");
      vnames.push_back("C_new");
      vecs.push_back(S_->GetPtr<CompositeVector>(conserved_key_, tag_current_).ptr());
      vecs.push_back(S_->GetPtr<CompositeVector>(conserved_key_, tag_next_).ptr());
      db_->WriteVectors(vnames, vecs, true);
    }
  } else {
    g->Update(1.0 / dt, *u_new, -1.0 / dt, *u_old, 0.0);
  }

  db_->WriteDivider();
  db_->WriteVector("res(acc)", g->Data().ptr());

  S_->GetEvaluator(cell_vol_key_, tag_next_).Update(*S_, name_);
  auto& cv = S_->Get<CompositeVector>(cell_vol_key_, tag_next_);

  if (is_source_) {
    if (theta_ < 1.0) {
      //S_->GetEvaluator(source_key_, tag_current_).Update(*S_, name_);
      g->Data()->Multiply(
        -(1.0 - theta_), S_->Get<CompositeVector>(source_key_, tag_current_), cv, 1.);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        db_->WriteVector(
          "source0", S_->GetPtr<CompositeVector>(source_key_, tag_current_).ptr(), false);
      }
    }

    if (theta_ > 0.0) {
      S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);
      g->Data()->Multiply(-theta_, S_->Get<CompositeVector>(source_key_, tag_next_), cv, 1.);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        db_->WriteVector(
          "source1", S_->GetPtr<CompositeVector>(source_key_, tag_next_).ptr(), false);
      }
    }
    db_->WriteVector("res(source)", g->Data().ptr());
  }
}


// updates the preconditioner
void
SurfaceBalanceBase::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_Physical_Default::Solution_to_State(*up, tag_next_);

  if (conserved_quantity_) {
    preconditioner_->Init();

    // add derivative of conserved quantity wrt primary
    S_->GetEvaluator(conserved_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
    auto dconserved_dT =
      S_->GetDerivativePtr<CompositeVector>(conserved_key_, tag_next_, key_, tag_next_);
    db_->WriteVector("d(cons)/d(prim)", dconserved_dT.ptr());
    preconditioner_acc_->AddAccumulationTerm(*dconserved_dT, h, "cell", false);

    // add derivative of source wrt primary
    if (theta_ > 0.0 && is_source_ && is_source_differentiable_ &&
        S_->GetEvaluator(source_key_, tag_next_).IsDifferentiableWRT(*S_, key_, tag_next_)) {
      Teuchos::RCP<const CompositeVector> dsource_dT;
      if (!source_finite_difference_) {
        // evaluate the derivative through chain rule and the DAG
        S_->GetEvaluator(source_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
        dsource_dT = S_->GetDerivativePtr<CompositeVector>(source_key_, tag_next_, key_, tag_next_);
      } else {
        // evaluate the derivative through finite differences
        S_->GetW<CompositeVector>(key_, tag_next_, name_).Shift(eps_);
        ChangedSolution();
        S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);
        auto dsource_dT_nc =
          Teuchos::rcp(new CompositeVector(S_->Get<CompositeVector>(source_key_, tag_next_)));

        S_->GetW<CompositeVector>(key_, tag_next_, name_).Shift(-eps_);
        ChangedSolution();
        S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);

        dsource_dT_nc->Update(
          -1 / eps_, S_->Get<CompositeVector>(source_key_, tag_next_), 1 / eps_);
        dsource_dT = dsource_dT_nc;
      }

      db_->WriteVector("d(Q)/d(prim)", dsource_dT.ptr());
      preconditioner_acc_->AddAccumulationTerm(*dsource_dT, -1.0 / theta_, "cell", true);
    }
  }
}


// applies preconditioner to u and returns the result in Pu
int
SurfaceBalanceBase::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                        Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon application:" << std::endl;

  if (conserved_quantity_) {
    db_->WriteVector("seb_res", u->Data().ptr(), true);
    preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
    db_->WriteVector("PC*p_res", Pu->Data().ptr(), true);
  } else {
    *Pu = *u;
    Pu->Scale(S_->get_time(tag_next_) - S_->get_time(tag_current_));
  }
  return 0;
}


bool
SurfaceBalanceBase::ModifyPredictor(double h,
                                    Teuchos::RCP<const TreeVector> u0,
                                    Teuchos::RCP<TreeVector> u)
{
  if (modify_predictor_positivity_preserving_) {
    int n_modified = 0;
    for (auto compname : *u->Data()) {
      auto& u_c = *u->Data()->ViewComponent(compname, false);
      for (int i = 0; i != u_c.MyLength(); ++i) {
        if (u_c[0][i] < 0.) {
          n_modified++;
          u_c[0][i] = 0.;
        }
      }
    }
    int n_modified_global = 0;
    u->Data()->Comm()->SumAll(&n_modified, &n_modified_global, 1);
    return n_modified_global > 0;
  }
  return false;
}

} // namespace SurfaceBalance
} // namespace Amanzi
