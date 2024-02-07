/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "CarbonSimple.hh"

namespace Amanzi {
namespace BGC {

CarbonSimple::CarbonSimple(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution)
  : Amanzi::PK(pk_tree, glist, S, solution),
    Amanzi::PK_Physical_Explicit_Default(pk_tree, glist, S, solution),
    is_diffusion_(false),
    is_source_(false),
    is_decomp_(false),
    npools_(-1)
{
  is_diffusion_ = plist_->get<bool>("is cryoturbation", true);
  if (is_diffusion_)
    div_diff_flux_key_ =
      Keys::readKey(*plist_, domain_, "divergence of bioturbation fluxes", "div_bioturbation");

  is_source_ = plist_->get<bool>("is source", true);
  if (is_source_) source_key_ = Keys::readKey(*plist_, domain_, "source", "carbon_source");

  is_decomp_ = plist_->get<bool>("is decomposition", true);
  if (is_decomp_)
    decomp_key_ =
      Keys::readKey(*plist_, domain_, "decomposition rate", "carbon_decomposition_rate");
}


// Setup data
void
CarbonSimple::Setup()
{
  PK_Physical_Explicit_Default::Setup();

  // number of carbon pools
  npools_ = plist_->get<int>("number of carbon pools");

  // cell volume
  if (cell_vol_key_ == std::string()) {
    cell_vol_key_ = plist_->get<std::string>("cell volume key", "cell_volume");
  }
  S_->Require<CompositeVector, CompositeVectorSpace>(cell_vol_key_, tag_current_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(cell_vol_key_, tag_current_);

  // diffusion
  if (is_diffusion_) {
    S_->Require<CompositeVector, CompositeVectorSpace>(div_diff_flux_key_, tag_current_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, npools_);
    S_->RequireEvaluator(div_diff_flux_key_, tag_current_);
  }

  // source terms
  if (is_source_) {
    S_->Require<CompositeVector, CompositeVectorSpace>(source_key_, tag_current_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, npools_);
    S_->RequireEvaluator(source_key_, tag_current_);
  }

  // decomposition terms
  if (is_decomp_) {
    S_->Require<CompositeVector, CompositeVectorSpace>(decomp_key_, tag_current_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, npools_);
    S_->RequireEvaluator(decomp_key_, tag_current_);
  }
}


// computes the non-linear functional f = f(t,u,udot)
void
CarbonSimple::FunctionalTimeDerivative(const double t, const TreeVector& u, TreeVector& f)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  // eventually we need to ditch this multi-state approach --etc
  AMANZI_ASSERT(std::abs(S_->get_time(tag_current_) - t) <
                1.e-4 * S_->get_time(tag_next_) - S_->get_time(tag_current_));
  PK_Physical_Default::Solution_to_State(u, tag_current_);

  // debugging
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Explicit deriv calculation: t = " << t << std::endl;
  db_->WriteCellInfo(true);
  db_->WriteVector("C_old", S_->GetPtr<CompositeVector>(key_, tag_current_).ptr());

  // Evaluate the derivative
  Teuchos::RCP<CompositeVector> dudt = f.Data();

  // -- apply the diffusion operator for cryoturbation
  ApplyDiffusion_(dudt.ptr());

  // -- add in source terms
  AddSources_(dudt.ptr());

  // -- add in decomposition
  AddDecomposition_(dudt.ptr());

  // scale all by cell volume
  S_->GetEvaluator(cell_vol_key_, tag_current_).Update(*S_, name_);
  const Epetra_MultiVector& cv =
    *S_->Get<CompositeVector>(cell_vol_key_, tag_current_).ViewComponent("cell", false);
  Epetra_MultiVector& dudt_c = *dudt->ViewComponent("cell", false);
  for (int c = 0; c != dudt_c.MyLength(); ++c) { dudt_c[0][c] *= cv[0][c]; }
}


// -- Calculate any diagnostics prior to doing vis
void
CarbonSimple::CalculateDiagnostics(const Tag& tag)
{
  // Call the functional.  This ensures that the vis gets updated values,
  // despite the fact that they have not yet been updated.
  TreeVector dudt(*solution_);
  FunctionalTimeDerivative(S_->get_time(tag_current_), *solution_old_, dudt);
}


// Physical routine to apply cryoturbation.
void
CarbonSimple::ApplyDiffusion_(const Teuchos::Ptr<CompositeVector>& g)
{
  if (is_diffusion_) {
    S_->GetEvaluator(div_diff_flux_key_, tag_current_).Update(*S_, name_);
    auto diff = S_->GetPtr<CompositeVector>(div_diff_flux_key_, tag_current_);
    g->Update(1., *diff, 0.);
    db_->WriteVector(" turbation rate", diff.ptr(), true);
  } else {
    g->PutScalar(0.);
  }
}

// Add in sources
void
CarbonSimple::AddSources_(const Teuchos::Ptr<CompositeVector>& g)
{
  if (is_source_) {
    S_->GetEvaluator(source_key_, tag_current_).Update(*S_, name_);
    auto src = S_->GetPtr<CompositeVector>(source_key_, tag_current_);
    g->Update(1., *src, 1.);
    db_->WriteVector(" source", src.ptr(), true);
  }
}

// Add in decomp
void
CarbonSimple::AddDecomposition_(const Teuchos::Ptr<CompositeVector>& g)
{
  if (is_decomp_) {
    S_->GetEvaluator(decomp_key_, tag_current_).Update(*S_, name_);
    auto src = S_->GetPtr<CompositeVector>(decomp_key_, tag_current_);
    g->Update(1., *src, 1.);
    db_->WriteVector(" decomp", src.ptr(), true);
  }
}


} // namespace BGC
} // namespace Amanzi
