/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "snow_distribution.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void
SnowDistribution::ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(tag);
  auto cond = S_->GetPtrW<CompositeVector>(uw_cond_key_, tag, name_);

  // update the stiffness matrix
  matrix_->Init();
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  matrix_diff_->ApplyBCs(true, true, true);

  // update the potential
  S_->GetEvaluator(potential_key_, tag).Update(*S_, name_);
  auto potential = S_->GetPtr<CompositeVector>(potential_key_, tag);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*potential, *g);
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void
SnowDistribution::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  // get these fields
  Teuchos::RCP<const CompositeVector> h1 = S_->GetPtr<CompositeVector>(key_, tag_next_);
  Epetra_MultiVector h1_positive(*h1->ViewComponent("cell", false));
  const auto& h1_v(*h1->ViewComponent("cell", false));
  for (int c = 0; c != h1_positive.MyLength(); ++c)
    h1_positive[0][c] = h1_v[0][c] > 0. ? h1_v[0][c] : 0.;

  Teuchos::RCP<const CompositeVector> h0 = S_->GetPtr<CompositeVector>(key_, tag_current_);
  Epetra_MultiVector h0_positive(*h0->ViewComponent("cell", false));
  const auto& h0_v(*h0->ViewComponent("cell", false));
  for (int c = 0; c != h0_positive.MyLength(); ++c)
    h0_positive[0][c] = h0_v[0][c] > 0. ? h0_v[0][c] : 0.;

  Teuchos::RCP<const CompositeVector> cv1 = S_->GetPtr<CompositeVector>(cv_key_, tag_next_);

  // note 10 is for conversion from precip m SWE to actual m
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);
  g->ViewComponent("cell", false)
    ->Multiply(10 * dt_factor_ / dt, *cv1->ViewComponent("cell", false), h1_positive, 1.);
  g->ViewComponent("cell", false)
    ->Multiply(-10 * dt_factor_ / dt, *cv1->ViewComponent("cell", false), h0_positive, 1.);
}

} // namespace Flow
} // namespace Amanzi
