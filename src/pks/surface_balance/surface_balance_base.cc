/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Adam Atchley, Satish Karra

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */
#include "LinearOperatorFactory.hh"
#include "surface_balance_base.hh"
#include "LinearOperatorFactory.hh"

namespace Amanzi {
namespace SurfaceBalance {


SurfaceBalanceBase::SurfaceBalanceBase(Teuchos::ParameterList& pk_tree,
                                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                       const Teuchos::RCP<State>& S,
                                       const Teuchos::RCP<TreeVector>& solution):
    PK(pk_tree, global_list,  S, solution),
    PK_PhysicalBDF_Default(pk_tree, global_list,  S, solution)
{
  // name the layer
  layer_ = plist_->get<std::string>("layer name", domain_);

  // source terms
  is_source_ = plist_->get<bool>("source term", true);
  if (is_source_) {
    source_key_ = Keys::readKey(*plist_, layer_, "source", "source_sink");
  }
  is_source_differentiable_ = plist_->get<bool>("source term is differentiable", true);
  source_finite_difference_ = plist_->get<bool>("source term finite difference", false);
  eps_ = plist_->get<double>("source term finite difference epsilon", 1.e-8);

  modify_predictor_positivity_preserving_ = plist_->get<bool>("modify predictor positivity preserving", false);
  
  theta_ = plist_->get<double>("time discretization theta", 1.0);
  AMANZI_ASSERT(theta_ <= 1.);
  AMANZI_ASSERT(theta_ >= 0.);

  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .01 * 55000.); // h * nl
}

// main methods
// -- Setup data.
void
SurfaceBalanceBase::Setup(const Teuchos::Ptr<State>& S) {
  PK_PhysicalBDF_Default::Setup(S);

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);

  // requirements: source terms from above
  if (is_source_) {
    S->RequireField(source_key_)->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(source_key_);
  }
  
  conserved_quantity_ = conserved_key_ != key_;
  if (conserved_quantity_) {
    S->RequireField(conserved_key_)->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(conserved_key_);
  }

  // operator for inverse
  Teuchos::ParameterList& acc_plist = plist_->sublist("accumulation preconditioner");
  acc_plist.set("entity kind", "cell");
  preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_plist, mesh_));
  preconditioner_ = preconditioner_acc_->global_operator();

  //    symbolic assemble
  precon_used_ = plist_->isSublist("preconditioner");
  if (precon_used_) {
    preconditioner_->SymbolicAssembleMatrix();
    preconditioner_->InitializePreconditioner(plist_->sublist("preconditioner"));
  }

  //    Potentially create a linear solver
  if (plist_->isSublist("linear solver")) {
    Teuchos::ParameterList& linsolve_sublist = plist_->sublist("linear solver");
    if (!linsolve_sublist.isSublist("verbose object"))
      linsolve_sublist.set("verbose object", plist_->sublist("verbose object"));
    AmanziSolvers::LinearOperatorFactory<Operators::Operator,CompositeVector,CompositeVectorSpace> fac;

    lin_solver_ = fac.Create(linsolve_sublist, preconditioner_);
  } else {
    lin_solver_ = preconditioner_;
  }
}


// computes the non-linear functional g = g(t,u,udot)
void
SurfaceBalanceBase::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                            Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  Teuchos::OSTab tab = vo_->getOSTab();
  double dt = t_new - t_old;

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, S_next_);
  
  bool debug = false;
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) debug = true;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << dt << std::endl;
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("u_old"); vnames.push_back("u_new");
    vecs.push_back(S_inter_->GetFieldData(key_).ptr());
    vecs.push_back(u_new->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  if (conserved_quantity_) {
    S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
    S_inter_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_inter_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> conserved1 = S_next_->GetFieldData(conserved_key_);
    Teuchos::RCP<const CompositeVector> conserved0 = S_inter_->GetFieldData(conserved_key_);
    g->Data()->Update(1.0/dt, *conserved1, -1.0/dt, *conserved0, 0.0);

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      std::vector<std::string> vnames;
      std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
      vnames.push_back("C_old"); vnames.push_back("C_new");
      vecs.push_back(S_inter_->GetFieldData(conserved_key_).ptr());
      vecs.push_back(S_next_->GetFieldData(conserved_key_).ptr());
      db_->WriteVectors(vnames, vecs, true);
    }
  } else {
    g->Update(1.0/dt, *u_new, -1.0/dt, *u_old, 0.0);
  }

  db_->WriteDivider();
  db_->WriteVector("res(acc)", g->Data().ptr());

  S_next_->GetFieldEvaluator(cell_vol_key_)->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> cv = S_next_->GetFieldData(cell_vol_key_);

  if (is_source_) {
    if (theta_ < 1.0) {
      S_inter_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_inter_.ptr(), name_);
      g->Data()->Multiply(-(1.0 - theta_), *S_inter_->GetFieldData(source_key_), *cv, 1.);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        db_->WriteVector("source0", S_inter_->GetFieldData(source_key_).ptr(), false);
      }
    }
    if (theta_ > 0.0) {
      S_next_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_next_.ptr(), name_);
      g->Data()->Multiply(-theta_, *S_next_->GetFieldData(source_key_), *cv, 1.);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        db_->WriteVector("source1", S_next_->GetFieldData(source_key_).ptr(), false);
      }
    }
    db_->WriteVector("res(source)", g->Data().ptr());
  }
}

  
// updates the preconditioner
void
SurfaceBalanceBase::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PK_Physical_Default::Solution_to_State(*up, S_next_);

  if (conserved_quantity_) {
    preconditioner_->Init();

    // add derivative of conserved quantity wrt primary
    S_next_->GetFieldEvaluator(conserved_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    auto dconserved_dT = S_next_->GetFieldData(Keys::getDerivKey(conserved_key_, key_));
    db_->WriteVector("d(cons)/d(prim)", dconserved_dT.ptr());
    preconditioner_acc_->AddAccumulationTerm(*dconserved_dT, h, "cell", false);

    // add derivative of source wrt primary
    if (theta_ > 0. && is_source_ && is_source_differentiable_ &&
        S_next_->GetFieldEvaluator(source_key_)->IsDependency(S_next_.ptr(), key_)) {

      Teuchos::RCP<const CompositeVector> dsource_dT;
      if (!source_finite_difference_) {
        // evaluate the derivative through chain rule and the DAG
        S_next_->GetFieldEvaluator(source_key_)->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
        dsource_dT = S_next_->GetFieldData(Keys::getDerivKey(source_key_,key_));
      } else {
        // evaluate the derivative through finite differences
        S_next_->GetFieldData(key_, name_)->Shift(eps_);
        ChangedSolution();
        S_next_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_next_.ptr(), name_);
        auto dsource_dT_nc = Teuchos::rcp(new CompositeVector(*S_next_->GetFieldData(source_key_)));

        S_next_->GetFieldData(key_, name_)->Shift(-eps_);
        ChangedSolution();
        S_next_->GetFieldEvaluator(source_key_)->HasFieldChanged(S_next_.ptr(), name_);

        dsource_dT_nc->Update(-1/eps_, *S_next_->GetFieldData(source_key_), 1/eps_);
        dsource_dT = dsource_dT_nc;
      }

      db_->WriteVector("d(Q)/d(prim)", dsource_dT.ptr());
      preconditioner_acc_->AddAccumulationTerm(*dsource_dT, -1.0/theta_, "cell", true);
    }

    if (precon_used_) {
      preconditioner_->AssembleMatrix();
      preconditioner_->UpdatePreconditioner();
    }
  }
}


// applies preconditioner to u and returns the result in Pu
int SurfaceBalanceBase::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

  if (conserved_quantity_) {
    db_->WriteVector("seb_res", u->Data().ptr(), true);
    lin_solver_->ApplyInverse(*u->Data(), *Pu->Data());
    db_->WriteVector("PC*p_res", Pu->Data().ptr(), true);
  } else {
    *Pu = *u;
    Pu->Scale(S_next_->time() - S_inter_->time());
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
      auto& u_c = *u->Data()->ViewComponent(compname,false);
      for (int i=0; i!=u_c.MyLength(); ++i) {
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

} // namespace
} // namespace
