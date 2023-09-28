/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "pk_helpers.hh"
#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void
OverlandPressureFlow::ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  auto& markers = bc_markers();
  auto& values = bc_values();
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(tag);
  auto cond = S_->GetPtrW<CompositeVector>(uw_cond_key_, tag, name_);

  // update the stiffness matrix
  matrix_->Init();
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  FixBCsForOperator_(tag, matrix_diff_.ptr()); // deals with zero gradient case
  matrix_diff_->ApplyBCs(true, true, true);

  // derive fluxes -- this gets done independently fo update as precon does
  // not calculate fluxes.
  S_->GetEvaluator(potential_key_, tag).Update(*S_, name_);
  auto pres_elev = S_->GetPtr<CompositeVector>(potential_key_, tag);
  auto flux = S_->GetPtrW<CompositeVector>(flux_key_, tag, name_);
  matrix_diff_->UpdateFlux(pres_elev.ptr(), flux.ptr());
  changedEvaluatorPrimary(flux_key_, tag, *S_);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, *g);
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void
OverlandPressureFlow::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  // get these fields
  S_->GetEvaluator(conserved_key_, tag_next_).Update(*S_, name_);
  //  S_->GetEvaluator(conserved_key_, tag_current_).Update(*S_, name_); // for the future...
  Teuchos::RCP<const CompositeVector> wc1 = S_->GetPtr<CompositeVector>(conserved_key_, tag_next_);
  Teuchos::RCP<const CompositeVector> wc0 =
    S_->GetPtr<CompositeVector>(conserved_key_, tag_current_);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
    vnames.push_back("  WC_old");
    vnames.push_back("  WC_new");
    vecs.push_back(wc0.ptr());
    vecs.push_back(wc1.ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  // Water content only has cells, while the residual has cells and faces.
  g->ViewComponent("cell", false)
    ->Update(1.0 / dt,
             *wc1->ViewComponent("cell", false),
             -1.0 / dt,
             *wc0->ViewComponent("cell", false),
             1.0);
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void
OverlandPressureFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g)
{
  Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

  S_->GetEvaluator(cv_key_, tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& cv1 =
    *S_->Get<CompositeVector>(cv_key_, tag_next_).ViewComponent("cell", false);

  if (is_source_term_) {
    // Add in external source term.
    S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& source1 =
      *S_->Get<CompositeVector>(source_key_, tag_next_).ViewComponent("cell", false);
    db_->WriteVector("  source", S_->GetPtr<CompositeVector>(source_key_, tag_next_).ptr(), false);
    int ncells = g_c.MyLength();
    for (int c = 0; c != ncells; ++c) { g_c[0][c] -= cv1[0][c] * source1[0][c]; }
  }

  if (coupled_to_subsurface_via_head_) {
    // Add in source term from coupling.
    S_->GetEvaluator(ss_flux_key_, tag_next_).Update(*S_, name_);
    Teuchos::RCP<const CompositeVector> source1 =
      S_->GetPtr<CompositeVector>(ss_flux_key_, tag_next_);

    // source term is in units of [mol / s]
    g_c.Update(-1., *source1->ViewComponent("cell", false), 1.);
  }
};


// -------------------------------------------------------------
// Nonlinear source terms contribute to PC
// -------------------------------------------------------------
void
OverlandPressureFlow::AddSourcesToPrecon_(double h)
{
  // -- update the source term derivatives
  if (is_source_term_ && source_term_is_differentiable_ &&
      S_->GetEvaluator(source_key_, tag_next_).IsDifferentiableWRT(*S_, key_, tag_next_)) {
    S_->GetEvaluator(source_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
    preconditioner_acc_->AddAccumulationTerm(
      S_->GetDerivative<CompositeVector>(source_key_, tag_next_, key_, tag_next_),
      -1.0,
      "cell",
      true);
  }
}


} // namespace Flow
} // namespace Amanzi
