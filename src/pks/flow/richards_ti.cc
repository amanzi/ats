/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Op.hh"
#include "TensorVector.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void
Richards::FunctionalResidual(double t_old,
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

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, tag_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old << " t1 = " << t_new << " h = " << h
               << std::endl;

  // debugging -- write primary variables to screen
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames{ "p_old", "p_new" };
  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.emplace_back(S_->GetPtr<CompositeVector>(key_, tag_current_).ptr());
  vecs.emplace_back(u.ptr());
  db_->WriteVectors(vnames, vecs, true);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);
  db_->WriteBoundaryConditions(bc_markers(), bc_values());

  // diffusion term, treated implicitly
  ApplyDiffusion_(tag_next_, res.ptr());
  // if (vapor_diffusion_) AddVaporDiffusionResidual_(tag_next_, res.ptr());

  // more debugging -- write diffusion/flux variables to screen
  vnames = { "sl_old", "sl_new" };
  vecs = { S_->GetPtr<CompositeVector>(sat_key_, tag_current_).ptr(),
           S_->GetPtr<CompositeVector>(sat_key_, tag_next_).ptr() };

  if (S_->HasRecordSet(sat_ice_key_)) {
    vnames.emplace_back("si_old");
    vnames.emplace_back("si_new");
    vecs.emplace_back(
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "saturation_ice"), tag_current_).ptr());
    vecs.emplace_back(
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "saturation_ice"), tag_next_).ptr());
  }
  vnames.emplace_back("k_rel");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(coef_key_, tag_next_).ptr());
  vnames.emplace_back("wind");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(flux_dir_key_, tag_next_).ptr());
  vnames.emplace_back("uw_k_rel");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_).ptr());
  vnames.emplace_back("flux");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(flux_key_, tag_next_).ptr());
  db_->WriteVectors(vnames, vecs, true);
  db_->WriteVector("res (diff)", res.ptr(), true);

  // accumulation term
  AddAccumulation_(res.ptr());
  db_->WriteVector("res (acc)", res.ptr(), true);

  // more debugging -- write accumulation variables to screen
  vnames = { "poro", "WC_old", "WC_new" };
  vecs = { S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "porosity"), tag_next_).ptr(),
           S_->GetPtr<CompositeVector>(conserved_key_, tag_current_).ptr(),
           S_->GetPtr<CompositeVector>(conserved_key_, tag_next_).ptr() };
  db_->WriteVectors(vnames, vecs);
  db_->WriteVector("res (acc)", res.ptr(), true);

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
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int
Richards::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon application:" << std::endl;

  // Apply the preconditioner
  db_->WriteVector("p_res", u->Data().ptr(), true);
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
  db_->WriteVector("PC*p_res", Pu->Data().ptr(), true);

  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void
Richards::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon update at t = " << t << std::endl;

  // Recreate mass matrices
  if ((!deform_key_.empty() &&
       S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " precon")) ||
      S_->GetEvaluator(perm_key_, tag_next_).Update(*S_, name_+" precon"))
    preconditioner_diff_->SetTensorCoefficient(Teuchos::rcpFromRef(S_->Get<TensorVector>(perm_key_, tag_next_).data));

  // update state with the solution up.
  if (std::abs(t - iter_counter_time_) / t > 1.e-4) {
    iter_ = 0;
    iter_counter_time_ = t;
  }
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);

  // update the rel perm according to the scheme of choice, also upwind derivatives of rel perm
  UpdatePermeabilityData_(tag_next_);
  if (jacobian_ && iter_ >= jacobian_lag_) UpdatePermeabilityDerivativeData_(tag_next_);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);

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

  // update mass matrix?
  if (S_->GetEvaluator(perm_key_, tag_next_).Update(*S_, name_+" precon"))
    preconditioner_diff_->SetTensorCoefficient(Teuchos::rcpFromRef(S_->Get<TensorVector>(perm_key_, tag_next_).data));

  // -- local matries, primary term
  preconditioner_->Init();
  preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  preconditioner_diff_->ApplyBCs(true, true, true);

  // -- local matries, Jacobian term
  if (jacobian_ && iter_ >= jacobian_lag_) {
    Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(flux_key_, tag_next_, name_);
    preconditioner_diff_->UpdateFlux(up->Data().ptr(), flux.ptr());
    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), up->Data().ptr());
  }

  // Update the preconditioner with accumulation terms.
  // -- update the accumulation derivatives
  S_->GetEvaluator(conserved_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> dwc_dp =
    S_->GetDerivativePtr<CompositeVector>(conserved_key_, tag_next_, key_, tag_next_);
  db_->WriteVector("    dwc_dp", dwc_dp.ptr());

  // -- update the cell-cell block
  preconditioner_acc_->AddAccumulationTerm(*dwc_dp, h, "cell", false);

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(h);

  // increment the iterator count
  iter_++;
};


} // namespace Flow
} // namespace Amanzi
