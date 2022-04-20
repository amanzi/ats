/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (h + elev)
// -------------------------------------------------------------
void OverlandPressureFlow::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g)
{
  auto& markers = bc_markers();
  auto& values = bc_values();
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(S_next_.ptr());
  auto cond = S_next_->GetFieldData(uw_cond_key_, name_);

  // update the stiffness matrix
  matrix_->Init();
  matrix_diff_->SetScalarCoefficient(cond, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  FixBCsForOperator_(S_next_.ptr(), matrix_diff_.ptr()); // deals with zero gradient case
  matrix_diff_->ApplyBCs(true, true, true);

  // derive fluxes -- this gets done independently fo update as precon does
  // not calculate fluxes.
  S->GetFieldEvaluator(potential_key_)->HasFieldChanged(S.ptr(), name_);
  auto pres_elev = S->GetFieldData(potential_key_);
  auto flux = S->GetFieldData(flux_key_, name_);
  matrix_diff_->UpdateFlux(pres_elev.ptr(), flux.ptr());
  if (S == S_next_.ptr()) flux_pvfe_->SetFieldAsChanged(S);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres_elev, *g);
};


// -------------------------------------------------------------
// Accumulation of water, dh/dt
// -------------------------------------------------------------
void OverlandPressureFlow::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  double dt = S_next_->time() - S_inter_->time();

  // get these fields
  S_next_->GetFieldEvaluator(conserved_key_)
      ->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(conserved_key_)
      ->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> wc1 =
      S_next_->GetFieldData(conserved_key_);
  Teuchos::RCP<const CompositeVector> wc0 =
      S_inter_->GetFieldData(conserved_key_);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("  WC_old"); vnames.push_back("  WC_new");
    vecs.push_back(wc0.ptr());
    vecs.push_back(wc1.ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  // Water content only has cells, while the residual has cells and faces.
  g->ViewComponent("cell",false)->Update(1.0/dt, *wc1->ViewComponent("cell",false),
          -1.0/dt, *wc0->ViewComponent("cell",false), 1.0);
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandPressureFlow::AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g)
{
  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

  const Epetra_MultiVector& cv1 =
    *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false);

  if (is_source_term_) {
    // Add in external source term.
    S_next_->GetFieldEvaluator(source_key_)
        ->HasFieldChanged(S_next_.ptr(), name_);
    const Epetra_MultiVector& source1 =
        *S_next_->GetFieldData(source_key_)->ViewComponent("cell",false);
    db_->WriteVector("  source", S_next_->GetFieldData(source_key_).ptr(), false);

    if (source_in_meters_) {
      // External source term is in [m water / s], not in [mols / s], so a
      // density is required.  This density should be upwinded.
      S_next_->GetFieldEvaluator(molar_dens_key_)
          ->HasFieldChanged(S_next_.ptr(), name_);
      S_next_->GetFieldEvaluator(source_molar_dens_key_)
          ->HasFieldChanged(S_next_.ptr(), name_);

      const Epetra_MultiVector& nliq1 = *S_next_->GetFieldData(molar_dens_key_)
          ->ViewComponent("cell",false);
      const Epetra_MultiVector& nliq1_s = *S_next_->GetFieldData(source_molar_dens_key_)
          ->ViewComponent("cell",false);

      int ncells = g_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        double s1 = source1[0][c] > 0. ? source1[0][c] * nliq1_s[0][c] : source1[0][c] * nliq1[0][c];
        g_c[0][c] -= cv1[0][c] * s1;
      }
    } else {
      int ncells = g_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        g_c[0][c] -= cv1[0][c] * source1[0][c];
      }
    }
  }

  if (coupled_to_subsurface_via_head_) {
    // Add in source term from coupling.
    S_next_->GetFieldEvaluator(ss_flux_key_)
        ->HasFieldChanged(S_next_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> source1 = S_next_->GetFieldData(ss_flux_key_);

    // source term is in units of [mol / s]
    g_c.Update(-1., *source1->ViewComponent("cell",false), 1.);
  }
};


} //namespace
} //namespace
