/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void
OverlandFlow::ApplyDiffusion_(const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& g)
{
  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(S_next_.ptr());

  // update the stiffness matrix
  matrix_->Init();
  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "upwind_overland_conductivity"), name_);
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update the potential
  S->GetEvaluator(Keys::getKey(domain_, "pres_elev"))->HasFieldChanged(S.ptr(), name_);

  // Patch up BCs for zero-gradient
  FixBCsForOperator_(S_next_.ptr());

  // derive fluxes -- this gets done independently fo update as precon does
  // not calculate fluxes.
  Teuchos::RCP<const CompositeVector> pres_elev =
    S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "pres_elev"));
  Teuchos::RCP<CompositeVector> flux = S->GetPtrW<CompositeVector>("surface-water_flux", name_);
  matrix_diff_->UpdateFlux(pres_elev.ptr(), flux.ptr());

  // assemble the stiffness matrix
  matrix_diff_->ApplyBCs(true, true, true);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, *g);
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void
OverlandFlow::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_inter_);

  // get these fields
  S_next_->GetEvaluator(Keys::getKey(domain_, "ponded_depth"))
    ->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetEvaluator(Keys::getKey(domain_, "ponded_depth"))
    ->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> wc1 =
    S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "ponded_depth"));
  Teuchos::RCP<const CompositeVector> wc0 =
    S_inter_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "ponded_depth"));
  Teuchos::RCP<const CompositeVector> cv =
    S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "cell_volume"));

  // Water content only has cells, while the residual has cells and faces.
  g->ViewComponent("cell", false)
    ->Multiply(1.0 / dt, *wc1->ViewComponent("cell", false), *cv->ViewComponent("cell", false), 1.);
  g->ViewComponent("cell", false)
    ->Multiply(
      -1.0 / dt, *wc0->ViewComponent("cell", false), *cv->ViewComponent("cell", false), 1.);
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void
OverlandFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g)
{
  Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

  const Epetra_MultiVector& cv1 =
    *S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "cell_volume"))
       ->ViewComponent("cell", false);

  if (is_source_term_) {
    // Add in external source term.
    S_next_->GetEvaluator(source_key_)->HasFieldChanged(S_next_.ptr(), name_);
    const Epetra_MultiVector& source1 =
      *S_next_->Get<CompositeVector>(source_key_).ViewComponent("cell", false);

    g_c.Multiply(-1., source1, cv1, 1.);
  }
};


} // namespace Flow
} // namespace Amanzi
