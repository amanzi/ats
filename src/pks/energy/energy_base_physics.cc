/* -*-  mode++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Solves:

de/dt + q dot grad h = div Ke grad T + S?
------------------------------------------------------------------------- */

#include "advection.hh"
#include "FieldEvaluator.hh"
#include "energy_base.hh"
#include "Op.hh"
#include "pk_helpers.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Accumulation of energy term de/dt
// -------------------------------------------------------------
void EnergyBase::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // update the energy at both the old and new times.
  S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // get the energy at each time
  Teuchos::RCP<const CompositeVector> e1 = S_next_->GetFieldData(conserved_key_);
  Teuchos::RCP<const CompositeVector> e0 = S_inter_->GetFieldData(conserved_key_);

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  g->ViewComponent("cell", false)
    ->Update(1.0/dt, *e1->ViewComponent("cell", false),
          -1.0/dt, *e0->ViewComponent("cell", false), 1.0);
};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, q dot grad h.
// -------------------------------------------------------------
void EnergyBase::AddAdvection_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g, bool negate) {

  // set up the operator

  // NOTE: fluxes are a MOLAR flux by choice of the flow pk, i.e.
  // [flux] =  mol/s

  // NOTE: this will be the eventual way to ensure it is up to date,
  // but there is no FieldEvaluator for darcy flux yet.  When there
  // is, we can take the evaluation out of Flow::commit_state(),
  // but for now we'll leave it there and assume it has been updated. --etc
  //  S->GetFieldEvaluator(flux_key_)->HasFieldChanged(S.ptr(), name_);
  ApplyDirichletBCsToEnthalpy_(S);

  // debugging
  db_->WriteBoundaryConditions(bc_adv_->bc_model(), bc_adv_->bc_value());
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_key_);
  Teuchos::RCP<const CompositeVector> enth = S->GetFieldData(enthalpy_key_);
  db_->WriteVectors({" adv flux", " enthalpy"}, {flux.ptr(), enth.ptr()}, true);

  matrix_adv_->global_operator()->Init();
  matrix_adv_->Setup(*flux);
  matrix_adv_->SetBCs(bc_adv_, bc_adv_);
  matrix_adv_->UpdateMatrices(flux.ptr());
  matrix_adv_->ApplyBCs(false, true, false);
  matrix_adv_->global_operator()->ComputeNegativeResidual(*enth, *g, false);
}

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void EnergyBase::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g) {
  // update the thermal conductivity
  UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);

  // update the stiffness matrix
  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // update the flux
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData(energy_flux_key_, name_);
  matrix_diff_->UpdateFlux(temp.ptr(), flux.ptr());

  // calculate the residual
  matrix_diff_->global_operator()->ComputeNegativeResidual(*temp, *g);
};


// ---------------------------------------------------------------------
// Add in energy source, which are accumulated by a single evaluator.
// ---------------------------------------------------------------------
void EnergyBase::AddSources_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // external sources of energy
  if (is_source_term_) {
    Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

    // Update the source term
    S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
    const Epetra_MultiVector& source1 =
        *S->GetFieldData(source_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& cv =
      *S->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

    // Add into residual
    unsigned int ncells = g_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      g_c[0][c] -= source1[0][c] * cv[0][c];
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "Adding external source term" << std::endl;
    db_->WriteVector("  Q_ext", S->GetFieldData(source_key_).ptr(), false);
    db_->WriteVector("res (src)", g, false);
  }
}


void EnergyBase::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
  // external sources of energy (temperature dependent source)
  if (is_source_term_ && is_source_term_differentiable_ &&
      S->GetFieldEvaluator(source_key_)->IsDependency(S, key_)) {

    Teuchos::RCP<const CompositeVector> dsource_dT;
    if (!is_source_term_finite_differentiable_) {
      // evaluate the derivative through the dag
      S->GetFieldEvaluator(source_key_)->HasFieldDerivativeChanged(S, name_, key_);
      dsource_dT = S->GetFieldData(Keys::getDerivKey(source_key_, key_));
    } else {
      // evaluate the derivative through finite differences
      double eps = 1.e-8;
      S->GetFieldData(key_, name_)->Shift(eps);
      ChangedSolution();
      S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
      auto dsource_dT_nc = Teuchos::rcp(new CompositeVector(*S->GetFieldData(source_key_)));

      S->GetFieldData(key_, name_)->Shift(-eps);
      ChangedSolution();
      S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);

      dsource_dT_nc->Update(-1/eps, *S->GetFieldData(source_key_), 1/eps);
      dsource_dT = dsource_dT_nc;
    }
    db_->WriteVector("  dQ_ext/dT", dsource_dT.ptr(), false);
    preconditioner_acc_->AddAccumulationTerm(*dsource_dT, -1.0, "cell", true);
  }
}

// -------------------------------------------------------------
// Plug enthalpy into the boundary faces manually.
// This will be removed once boundary faces exist.
// -------------------------------------------------------------
void EnergyBase::ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S) {

  // in the diffusive flux condition, first update the boundary face temperatures for FV
  auto& T_vec = *S->GetFieldData(key_, name_);
  if (T_vec.HasComponent("boundary_face")) {
    Epetra_MultiVector& T_bf = *T_vec.ViewComponent("boundary_face", false);
    const Epetra_MultiVector& T_c = *T_vec.ViewComponent("cell", false);

    for (int bf=0; bf!=T_bf.MyLength(); ++bf) {
      AmanziMesh::Entity_ID f = getBoundaryFaceFace(*mesh_, bf);

      // NOTE: this should get refactored into a helper class, much like predictor_delegate_bc_flux
      // as this would be necessary to deal with general discretizations.  Note that this is not
      // needed in cases where boundary faces are already up to date (e.g. MFD, maybe NLFV?)
      if (bc_markers()[f] == Operators::OPERATOR_BC_NEUMANN &&
          bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
        // diffusive flux BC
        AmanziMesh::Entity_ID c = getFaceOnBoundaryInternalCell(*mesh_, f);
        const auto& Acc = matrix_diff_->local_op()->matrices_shadow[f];
        T_bf[0][bf] = (Acc(0,0)*T_c[0][c] - bc_values()[f]*mesh_->face_area(f)) / Acc(0,0);
      }
    }
  }

  // then put the boundary fluxes in faces for Dirichlet BCs.
  S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S, name_);

  const Epetra_MultiVector& enth_bf =
    *S->GetFieldData(enthalpy_key_)->ViewComponent("boundary_face",false);

  int nbfaces = enth_bf.MyLength();
  for (int bf=0; bf!=nbfaces; ++bf) {
    AmanziMesh::Entity_ID f = getBoundaryFaceFace(*mesh_, bf);

    if (bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
      bc_adv_->bc_value()[f] = enth_bf[0][bf];
    }
  }
}



} //namespace Energy
} //namespace Amanzi
