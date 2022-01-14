/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#include "CarbonSimple.hh"

namespace Amanzi {
namespace BGC {


CarbonSimple::CarbonSimple(Teuchos::ParameterList& pk_tree,
                            const Teuchos::RCP<Teuchos::ParameterList>& glist,
                            const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<TreeVector>& solution):
    Amanzi::PK(pk_tree, glist, S, solution),
    Amanzi::PK_Physical_Explicit_Default(pk_tree, glist, S, solution),
    is_diffusion_(false),
    is_source_(false),
    is_decomp_(false),
    npools_(-1)
{}


// Setup data
void
CarbonSimple::Setup(const Teuchos::Ptr<State>& S) {
  PK_Physical_Explicit_Default::Setup(S);

  // number of carbon pools
  npools_ = plist_->get<int>("number of carbon pools");
  
  // cell volume
  if (cell_vol_key_ == std::string()) {
    cell_vol_key_ = plist_->get<std::string>("cell volume key", "cell_volume");
  }
  S->Require<CompositeVector,CompositeVectorSpace>(cell_vol_key_, Tags::NEXT).SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(cell_vol_key_);
  
  // diffusion
  is_diffusion_ = plist_->get<bool>("is cryoturbation", true);
  if (is_diffusion_) {
    div_diff_flux_key_ = plist_->get<std::string>("divergence of bioturbation fluxes", "div_bioturbation");

    S->Require<CompositeVector,CompositeVectorSpace>(div_diff_flux_key_, Tags::NEXT).SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, npools_);
    S->RequireEvaluator(div_diff_flux_key_);
  }

  // source terms
  is_source_ = plist_->get<bool>("is source", true);
  if (is_source_) {
    source_key_ = plist_->get<std::string>("source key", "carbon_source");

    S->Require<CompositeVector,CompositeVectorSpace>(source_key_, Tags::NEXT).SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, npools_);
    S->RequireEvaluator(source_key_);
  }

  // decomposition terms
  is_decomp_ = plist_->get<bool>("is decomposition", true);
  if (is_decomp_) {
    decomp_key_ = plist_->get<std::string>("decomposition rate", "carbon_decomposition_rate");

    S->Require<CompositeVector,CompositeVectorSpace>(div_diff_flux_key_, Tags::NEXT).SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, npools_);
    S->RequireEvaluator(div_diff_flux_key_);
  }
}


// computes the non-linear functional f = f(t,u,udot)
void
CarbonSimple::FunctionalTimeDerivative(const double t, const TreeVector& u, TreeVector& f) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  // eventually we need to ditch this multi-state approach --etc
  AMANZI_ASSERT(std::abs(S_->get_time(tag_inter_) - t) < 1.e-4*S_->get_time(tag_next_) - S_->get_time(tag_inter_));
  PK_Physical_Default::Solution_to_State(u, S_inter_);

  // debugging
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Explicit deriv calculation: t = " << t << std::endl;
  db_->WriteCellInfo(true);
  db_->WriteVector("C_old", S_inter_->GetPtr<CompositeVector>(key_).ptr());

  // Evaluate the derivative
  Teuchos::RCP<CompositeVector> dudt = f.Data();

  // -- apply the diffusion operator for cryoturbation
  ApplyDiffusion_(S_inter_.ptr(), dudt.ptr());

  // -- add in source terms
  AddSources_(S_inter_.ptr(), dudt.ptr());

  // -- add in decomposition
  AddDecomposition_(S_inter_.ptr(), dudt.ptr());

  // scale all by cell volume
  const Epetra_MultiVector& cv = *S_inter_->GetPtr<CompositeVector>(cell_vol_key_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& dudt_c = *dudt->ViewComponent("cell",false);
  for (int c=0; c!=dudt_c.MyLength(); ++c) {
    dudt_c[0][c] *= cv[0][c];
  }
}


// -- Calculate any diagnostics prior to doing vis
void
CarbonSimple::CalculateDiagnostics(const Teuchos::RCP<State>& S) {
  // Call the functional.  This ensures that the vis gets updated values,
  // despite the fact that they have not yet been updated.
  TreeVector dudt(*solution_);
  FunctionalTimeDerivative(S->get_time(), *solution_old_, dudt);
}


// Physical routine to apply cryoturbation.
void
CarbonSimple::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  if (is_diffusion_) {
    S->GetEvaluator(div_diff_flux_key_)->HasFieldChanged(S, name_);
    Teuchos::RCP<const CompositeVector> diff = S->GetPtr<CompositeVector>(div_diff_flux_key_);
    g->Update(1., *diff, 0.);
    db_->WriteVector(" turbation rate", diff.ptr(), true);
  } else {
    g->PutScalar(0.);
  }
}

// Add in sources
void
CarbonSimple::AddSources_(const Teuchos::Ptr<State>& S,
                          const Teuchos::Ptr<CompositeVector>& g) {
  if (is_source_) {
    S->GetEvaluator(source_key_)->HasFieldChanged(S, name_);
    Teuchos::RCP<const CompositeVector> src = S->GetPtr<CompositeVector>(source_key_);
    g->Update(1., *src, 1.);
    db_->WriteVector(" source", src.ptr(), true);
  }
}

// Add in decomp
void
CarbonSimple::AddDecomposition_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  if (is_decomp_) {
    S->GetEvaluator(decomp_key_)->HasFieldChanged(S, name_);
    Teuchos::RCP<const CompositeVector> src = S->GetPtr<CompositeVector>(decomp_key_);
    g->Update(1., *src, 1.);
    db_->WriteVector(" decomp", src.ptr(), true);
  }
}




} // namespace BGC
} // namespace ATS
