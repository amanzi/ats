/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"

#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void
OverlandPressureFlow::ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(tag);

  // update the stiffness matrix
  matrix_->Zero();
  matrix_diff_->SetScalarCoefficient(S_->GetPtr<CompositeVector>(uw_cond_key_, tag), Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  FixBCsForOperator_(tag, matrix_diff_.ptr()); // deals with zero gradient case

  // derive fluxes -- this gets done independently fo update as precon does
  // not calculate fluxes.
  S_->GetEvaluator(potential_key_, tag).Update(*S_, name_);
  auto pres_elev = S_->GetPtr<CompositeVector>(potential_key_, tag);
  auto flux = S_->GetPtrW<CompositeVector>(flux_key_, tag, name_);
  matrix_diff_->UpdateFlux(pres_elev.ptr(), flux.ptr());
  PKHelpers::changedEvaluatorPrimary(flux_key_, tag, *S_);

  // note tpetra has swapped order -- ApplyBCs after flux...
  matrix_diff_->ApplyBCs(true, true, true);

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
  g->getComponent("cell", false)
    ->update(1.0 / dt,
             *wc1->getComponent("cell", false),
             -1.0 / dt,
             *wc0->getComponent("cell", false),
             1.0);
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void
OverlandPressureFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g)
{
  S_->GetEvaluator(cell_vol_key_, tag_next_).Update(*S_, name_);

  if (is_source_term_) {
    // Add in external source term.
    S_->GetEvaluator(source_key_, tag_next_).Update(*S_, name_);
    db_->WriteVector("  source", S_->GetPtr<CompositeVector>(source_key_, tag_next_).ptr(), false);

    g->getComponent("cell", false)
      ->elementWiseMultiply(
        -1.0,
        *S_->Get<CompositeVector>(cell_vol_key_, tag_next_)
           .getComponent("cell", false)
           ->getVector(0),
        *S_->Get<CompositeVector>(source_key_, tag_next_).getComponent("cell", false),
        1.0);
  }

  if (coupled_to_subsurface_via_head_) {
    // Add in source term from coupling.
    S_->GetEvaluator(ss_flux_key_, tag_next_).Update(*S_, name_);
    const CompositeVector& source1 = S_->Get<CompositeVector>(ss_flux_key_, tag_next_);

    // source term is in units of [mol / s]
    g->getComponent("cell", false)->update(-1., *source1.getComponent("cell", false), 1.);
  }
};


} // namespace Flow
} // namespace Amanzi
