/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

/*
Steady state solution of Richards equation

*/

#include "TensorVector.hh"
#include "richards_steadystate.hh"

namespace Amanzi {
namespace Flow {

RichardsSteadyState::RichardsSteadyState(Teuchos::ParameterList& pk_tree,
                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                         const Teuchos::RCP<State>& S,
                                         const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, glist, S, solution), Richards(pk_tree, glist, S, solution)
{}


// RichardsSteadyState is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void
RichardsSteadyState::FunctionalResidual(double t_old,
                                        double t_new,
                                        Teuchos::RCP<const TreeVector> u_old,
                                        Teuchos::RCP<TreeVector> u_new,
                                        Teuchos::RCP<TreeVector> g)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  double h = t_new - t_old;
  AMANZI_ASSERT(std::abs(S_->get_time(tag_current_) - t_old) < 1.e-4 * h);
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4 * h);

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, tag_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old << " t1 = " << t_new << " h = " << h
               << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("p_old");
  vnames.push_back("p_new");
  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.push_back(S_->GetPtr<CompositeVector>(key_, tag_current_).ptr());
  vecs.push_back(u.ptr());
  db_->WriteVectors(vnames, vecs, true);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);
  db_->WriteBoundaryConditions(bc_markers(), bc_values());

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(tag_next_, res.ptr());

  // evaulate water content, because otherwise it is never done.
  S_->GetEvaluator(conserved_key_, tag_next_).Update(*S_, name_);

  // dump s_old, s_new
  vnames[0] = "sl_old";
  vnames[1] = "sl_new";
  vecs[0] = S_->GetPtr<CompositeVector>(sat_key_, tag_current_).ptr();
  vecs[1] = S_->GetPtr<CompositeVector>(sat_key_, tag_next_).ptr();

  if (S_->HasRecordSet(sat_ice_key_)) {
    vnames.push_back("si_old");
    vnames.push_back("si_new");
    vecs.push_back(
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "saturation_ice"), tag_current_).ptr());
    vecs.push_back(
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "saturation_ice"), tag_next_).ptr());
  }
  vnames.push_back("poro");
  vecs.push_back(S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "porosity"), tag_next_).ptr());
  vnames.push_back("k_rel");
  vecs.push_back(S_->GetPtr<CompositeVector>(coef_key_, tag_next_).ptr());
  vnames.push_back("wind");
  vecs.push_back(S_->GetPtr<CompositeVector>(flux_dir_key_, tag_next_).ptr());
  vnames.push_back("uw_k_rel");
  vecs.push_back(S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_).ptr());
  vnames.push_back("flux");
  vecs.push_back(S_->GetPtr<CompositeVector>(flux_key_, tag_next_).ptr());
  db_->WriteVectors(vnames, vecs, true);

  db_->WriteVector("res (diff)", res.ptr(), true);

  // source term
  if (is_source_term_) {
    if (explicit_source_) {
      AddSources_(tag_current_, res.ptr());
    } else {
      AddSources_(tag_next_, res.ptr());
    }
    db_->WriteVector("res (src)", res.ptr(), false);
  }
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void
RichardsSteadyState::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Precon update at t = " << t << std::endl;
  }

  // Recreate mass matrices
  if ((!deform_key_.empty() &&
       S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " precon")) ||
      S_->GetEvaluator(perm_key_, tag_next_).Update(*S_, name_ + " precon"))
    preconditioner_diff_->SetTensorCoefficient(
      Teuchos::rcpFromRef(S_->Get<TensorVector>(perm_key_, tag_next_).data));

  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);

  // update the rel perm according to the scheme of choice, also upwind derivatives of rel perm
  UpdatePermeabilityData_(tag_next_);
  if (jacobian_ && iter_ >= jacobian_lag_) UpdatePermeabilityDerivativeData_(tag_next_);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);

  // Zero out the preconditioner and local matrices
  preconditioner_->Init();

  // fill local matrices
  // -- gravity fluxes
  S_->GetEvaluator(mass_dens_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag_next_);
  preconditioner_diff_->SetDensity(rho);

  // -- jacobian term
  Teuchos::RCP<const CompositeVector> dkrdp = Teuchos::null;
  if (jacobian_ && iter_ >= jacobian_lag_) {
    if (!duw_coef_key_.empty()) {
      dkrdp = S_->GetPtr<CompositeVector>(duw_coef_key_, tag_next_);
    } else {
      dkrdp = S_->GetDerivativePtr<CompositeVector>(coef_key_, tag_next_, key_, tag_next_);
    }
  }

  // -- primary term
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_);
  preconditioner_diff_->SetScalarCoefficient(rel_perm, dkrdp);

  // -- fill local matrices
  preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  preconditioner_diff_->ApplyBCs(true, true, true);

  // -- update with Jacobian terms
  if (jacobian_ && iter_ >= jacobian_lag_) {
    Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(flux_key_, tag_next_, name_);
    preconditioner_diff_->UpdateFlux(up->Data().ptr(), flux.ptr());
    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), up->Data().ptr());
  }

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(h);

  // increment the iterator count
  iter_++;
};


} // namespace Flow
} // namespace Amanzi
