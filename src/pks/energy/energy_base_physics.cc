/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Solves:

de/dt + q dot grad h = div Ke grad T + S?
------------------------------------------------------------------------- */

#include "advection.hh"
#include "Evaluator.hh"
#include "energy_base.hh"
#include "Op.hh"
#include "pk_helpers.hh"
#include "Mesh_Algorithms.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Accumulation of energy term de/dt
// -------------------------------------------------------------
void
EnergyBase::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  // update the energy at both the old and new times.
  S_->GetEvaluator(conserved_key_, tag_next_).Update(*S_, name_);
  // S_->GetEvaluator(conserved_key_, tag_current_).Update(*S_, name_); // for the future...

  // get the energy at each time
  Teuchos::RCP<const CompositeVector> e1 = S_->GetPtr<CompositeVector>(conserved_key_, tag_next_);
  Teuchos::RCP<const CompositeVector> e0 =
    S_->GetPtr<CompositeVector>(conserved_key_, tag_current_);

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  g->ViewComponent("cell", false)
    ->Update(1.0 / dt,
             *e1->ViewComponent("cell", false),
             -1.0 / dt,
             *e0->ViewComponent("cell", false),
             1.0);
};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, q dot grad h.
// -------------------------------------------------------------
void
EnergyBase::AddAdvection_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g, bool negate)
{
  // set up the operator

  // NOTE: fluxes are a MOLAR flux by choice of the flow pk, i.e.
  // [flux] =  mol/s

  // NOTE: this will be the eventual way to ensure it is up to date,
  // but there is no Evaluator for darcy flux yet.  When there
  // is, we can take the evaluation out of Flow::commit_state(),
  // but for now we'll leave it there and assume it has been updated. --etc
  //  S->GetEvaluator(flux_key_)->HasFieldChanged(S.ptr(), name_);
  ApplyDirichletBCsToEnthalpy_(tag);

  // debugging
  db_->WriteBoundaryConditions(bc_adv_->bc_model(), bc_adv_->bc_value());
  Teuchos::RCP<const CompositeVector> flux = S_->GetPtr<CompositeVector>(flux_key_, tag);
  Teuchos::RCP<const CompositeVector> enth = S_->GetPtr<CompositeVector>(enthalpy_key_, tag);
  db_->WriteVectors({ " adv flux", " enthalpy" }, { flux.ptr(), enth.ptr() }, true);

  matrix_adv_->global_operator()->Init();
  matrix_adv_->Setup(*flux);
  matrix_adv_->SetBCs(bc_adv_, bc_adv_);
  matrix_adv_->UpdateMatrices(flux.ptr());
  matrix_adv_->ApplyBCs(false, true, false);

  // update the flux
  Teuchos::RCP<CompositeVector> adv_energy =
    S_->GetPtrW<CompositeVector>(adv_energy_flux_key_, tag, name_);
  matrix_adv_->UpdateFlux(enth.ptr(), flux.ptr(), bc_adv_, adv_energy.ptr());
  changedEvaluatorPrimary(adv_energy_flux_key_, tag, *S_);

  matrix_adv_->global_operator()->ComputeNegativeResidual(*enth, *g, false);
}

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void
EnergyBase::ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  // force mass matrices to change
  if (!deform_key_.empty() &&
      S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " matrix"))
    matrix_diff_->SetTensorCoefficient(Teuchos::null);

  // update the thermal conductivity
  UpdateConductivityData_(tag);
  auto cond = S_->GetPtrW<CompositeVector>(uw_conductivity_key_, tag, name_);

  Teuchos::RCP<const CompositeVector> temp = S_->GetPtrW<CompositeVector>(key_, tag, name_);

  // update the stiffness matrix
  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // update the flux
  Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(energy_flux_key_, tag, name_);
  matrix_diff_->UpdateFlux(temp.ptr(), flux.ptr());
  changedEvaluatorPrimary(energy_flux_key_, tag, *S_);

  // calculate the residual
  matrix_diff_->global_operator()->ComputeNegativeResidual(*temp, *g);
};


// ---------------------------------------------------------------------
// Add in energy source, which are accumulated by a single evaluator.
// ---------------------------------------------------------------------
void
EnergyBase::AddSources_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

  S_->GetEvaluator(cell_vol_key_, tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& cv =
    *S_->Get<CompositeVector>(cell_vol_key_, tag_next_).ViewComponent("cell", false);

  // external sources of energy
  if (is_source_term_) {
    // Update the source term
    S_->GetEvaluator(source_key_, tag).Update(*S_, name_);
    const Epetra_MultiVector& source1 =
      *S_->Get<CompositeVector>(source_key_, tag).ViewComponent("cell", false);

    // Add into residual
    unsigned int ncells = g_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) { g_c[0][c] -= source1[0][c] * cv[0][c]; }

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Adding external source term" << std::endl;
    db_->WriteVector("  Q_ext", S_->GetPtr<CompositeVector>(source_key_, tag).ptr(), false);
    db_->WriteVector("res (src)", g, false);
  }
}


void
EnergyBase::AddSourcesToPrecon_(double h)
{
  // external sources of energy (temperature dependent source)
  if (is_source_term_ && is_source_term_differentiable_ &&
      S_->GetEvaluator(source_key_, tag_next_).IsDifferentiableWRT(*S_, key_, tag_next_)) {
    Teuchos::RCP<CompositeVector> dsource_dT;

    if (is_source_term_finite_differentiable_) {
      // evaluate the derivative through finite differences
      double eps = 1.e-8;
      S_->GetW<CompositeVector>(key_, tag_next_, name_).Shift(eps);
      ChangedSolution();
      S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);
      auto dsource_dT_nc =
        Teuchos::rcp(new CompositeVector(S_->Get<CompositeVector>(source_key_, tag_next_)));

      S_->GetW<CompositeVector>(key_, tag_next_, name_).Shift(-eps);
      ChangedSolution();
      S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);

      dsource_dT_nc->Update(-1 / eps, S_->Get<CompositeVector>(source_key_, tag_next_), 1 / eps);
      dsource_dT = dsource_dT_nc;

    } else {
      // evaluate the derivative through the dag
      S_->GetEvaluator(source_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
      dsource_dT = S_->GetDerivativePtrW<CompositeVector>(
        source_key_, tag_next_, key_, tag_next_, source_key_);
    }
    db_->WriteVector("  dQ_ext/dT", dsource_dT.ptr(), false);
    preconditioner_acc_->AddAccumulationTerm(*dsource_dT, -1.0, "cell", true);
  }
}

// -------------------------------------------------------------
// Plug enthalpy into the boundary faces manually.
// This will be removed once boundary faces exist.
// -------------------------------------------------------------
void
EnergyBase::ApplyDirichletBCsToEnthalpy_(const Tag& tag)
{
  // in the diffusive flux condition, first update the boundary face temperatures for FV
  auto& T_vec = *S_->GetPtrW<CompositeVector>(key_, tag, name_);
  if (T_vec.HasComponent("boundary_face")) {
    Epetra_MultiVector& T_bf = *T_vec.ViewComponent("boundary_face", false);
    const Epetra_MultiVector& T_c = *T_vec.ViewComponent("cell", false);

    for (int bf = 0; bf != T_bf.MyLength(); ++bf) {
      AmanziMesh::Entity_ID f = getBoundaryFaceFace(*mesh_, bf);

      // NOTE: this should get refactored into a helper class, much like predictor_delegate_bc_flux
      // as this would be necessary to deal with general discretizations.  Note that this is not
      // needed in cases where boundary faces are already up to date (e.g. MFD, maybe NLFV?)
      if (bc_markers()[f] == Operators::OPERATOR_BC_NEUMANN &&
          bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
        // diffusive flux BC
        AmanziMesh::Entity_ID c = getFaceOnBoundaryInternalCell(*mesh_, f);
        const auto& Acc = matrix_diff_->local_op()->matrices_shadow[f];
        T_bf[0][bf] = (Acc(0, 0) * T_c[0][c] - bc_values()[f] * mesh_->getFaceArea(f)) / Acc(0, 0);
      }
    }
  }

  // then put the boundary fluxes in faces for Dirichlet BCs.
  S_->GetEvaluator(enthalpy_key_, tag).Update(*S_, name_);

  const Epetra_MultiVector& enth_bf =
    *S_->Get<CompositeVector>(enthalpy_key_, tag).ViewComponent("boundary_face", false);

  int nbfaces = enth_bf.MyLength();
  for (int bf = 0; bf != nbfaces; ++bf) {
    AmanziMesh::Entity_ID f = getBoundaryFaceFace(*mesh_, bf);

    if (bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
      bc_adv_->bc_value()[f] = enth_bf[0][bf];
    }
  }
}


} //namespace Energy
} //namespace Amanzi
