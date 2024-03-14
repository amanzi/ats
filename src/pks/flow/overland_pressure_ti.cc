/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_pressure.hh"
#include "Op.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_RES_FLAG 0

// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void
OverlandPressureFlow::FunctionalResidual(double t_old,
                                         double t_new,
                                         Teuchos::RCP<TreeVector> u_old,
                                         Teuchos::RCP<TreeVector> u_new,
                                         Teuchos::RCP<TreeVector> g)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  // bookkeeping
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

  // unnecessary here if not debeugging, but doesn't hurt either
  S_->GetEvaluator(potential_key_, tag_next_).Update(*S_, name_);

  // debugging -- write primary variables to screen
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames{ "p_old", "p_new", "z", "h_old", "h_new", "h+z" };
  if (plist_->isSublist("overland conductivity subgrid evaluator")) {
    vnames.emplace_back("pd - dd");
    vnames.emplace_back("frac_cond");
  }
  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.emplace_back(S_->GetPtr<CompositeVector>(key_, tag_current_).ptr());
  vecs.emplace_back(u.ptr());
  vecs.emplace_back(S_->GetPtr<CompositeVector>(elev_key_, tag_next_).ptr());
  vecs.emplace_back(S_->GetPtr<CompositeVector>(pd_key_, tag_current_).ptr());
  vecs.emplace_back(S_->GetPtr<CompositeVector>(pd_key_, tag_next_).ptr());
  vecs.emplace_back(S_->GetPtr<CompositeVector>(potential_key_, tag_next_).ptr());
  if (plist_->isSublist("overland conductivity subgrid evaluator")) {
    // fixme -- add keys!
    vecs.emplace_back(
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "mobile_depth"), tag_next_).ptr());
    vecs.emplace_back(
      S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "fractional_conductance"), tag_next_)
        .ptr());
  }
  db_->WriteVectors(vnames, vecs, true);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);
  db_->WriteBoundaryConditions(bc_markers(), bc_values());

  // diffusion term, treated implicitly
  ApplyDiffusion_(tag_next_, res.ptr());

  // more debugging -- write diffusion/flux variables to screen
  vnames.clear();
  vecs.clear();
  if (S_->HasRecord(Keys::getKey(domain_, "unfrozen_fraction"), tag_next_) &&
      S_->HasRecord(Keys::getKey(domain_, "unfrozen_fraction"), tag_current_)) {
    Key uf_key = Keys::getKey(domain_, "unfrozen_fraction");
    vnames = { "uf_frac_old", "uf_frac_new" };
    vecs = { S_->GetPtr<CompositeVector>(uf_key, tag_current_).ptr(),
             S_->GetPtr<CompositeVector>(uf_key, tag_next_).ptr() };
  }
  vnames.emplace_back("uw_dir");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(flux_dir_key_, tag_next_).ptr());
  vnames.emplace_back("k");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(cond_key_, tag_next_).ptr());
  vnames.emplace_back("k_uw");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(uw_cond_key_, tag_next_).ptr());
  vnames.emplace_back("q_surf");
  vecs.emplace_back(S_->GetPtr<CompositeVector>(flux_key_, tag_next_).ptr());
  db_->WriteVectors(vnames, vecs, true);
  db_->WriteVector("res (diff)", res.ptr(), true);

  // accumulation term
  AddAccumulation_(res.ptr());
  db_->WriteVector("res (acc)", res.ptr(), true);

  // add rhs load value
  AddSourceTerms_(res.ptr());
  db_->WriteVector("res (src)", res.ptr(), true);
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int
OverlandPressureFlow::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                          Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon application:" << std::endl;
  AMANZI_ASSERT(!precon_scaled_); // otherwise this factor was built into the matrix

  // apply the preconditioner
  db_->WriteVector("h_res", u->Data().ptr(), true);
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
  db_->WriteVector("PC*h_res (h-coords)", Pu->Data().ptr(), true);

  // tack on the variable change
  const Epetra_MultiVector& dh_dp =
    *S_->GetDerivative<CompositeVector>(pd_bar_key_, tag_next_, key_, tag_next_)
       .ViewComponent("cell", false);
  Epetra_MultiVector& Pu_c = *Pu->Data()->ViewComponent("cell", false);

  unsigned int ncells = Pu_c.MyLength();
  for (unsigned int c = 0; c != ncells; ++c) { Pu_c[0][c] /= dh_dp[0][c]; }
  db_->WriteVector("PC*h_res (p-coords)", Pu->Data().ptr(), true);

  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void
OverlandPressureFlow::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  if (std::abs(t - iter_counter_time_) / t > 1.e-4) {
    iter_ = 0;
    iter_counter_time_ = t;
  }

  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);

  // calculating the operator is done in 3 steps:
  // 1. Diffusion components

  // 1.a: Pre-assembly updates.
  // -- update the rel perm according to the boundary info and upwinding
  // -- scheme of choice
  UpdatePermeabilityData_(tag_next_);
  if (jacobian_ && iter_ >= jacobian_lag_) UpdatePermeabilityDerivativeData_(tag_next_);

  // -- update boundary condition markers, which set the BC type
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);

  // fill local matrices
  // -- jacobian term
  Teuchos::RCP<const CompositeVector> dcond = Teuchos::null;
  if (jacobian_ && iter_ >= jacobian_lag_) {
    if (!duw_cond_key_.empty()) {
      dcond = S_->GetPtr<CompositeVector>(duw_cond_key_, tag_next_);
    } else {
      dcond = S_->GetDerivativePtr<CompositeVector>(cond_key_, tag_next_, pd_key_, tag_next_);
    }
  }
  // -- primary term
  Teuchos::RCP<const CompositeVector> cond = S_->GetPtr<CompositeVector>(uw_cond_key_, tag_next_);
  preconditioner_diff_->SetScalarCoefficient(cond, dcond);

  // -- local matrices, primary term
  preconditioner_->Init();
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // -- local matrices, Jacobian term
  if (jacobian_ && iter_ >= jacobian_lag_) {
    S_->GetEvaluator(potential_key_, tag_next_).Update(*S_, name_);
    Teuchos::RCP<const CompositeVector> pres_elev =
      S_->GetPtr<CompositeVector>(potential_key_, tag_next_);
    Teuchos::RCP<CompositeVector> flux = Teuchos::null;
    if (preconditioner_->RangeMap().HasComponent("face")) {
      flux = S_->GetPtrW<CompositeVector>(flux_key_, tag_next_, name_);
      preconditioner_diff_->UpdateFlux(pres_elev.ptr(), flux.ptr());
    }
    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), pres_elev.ptr());
  }

  // 2. Accumulation shift
  //    The desire is to keep this matrix invertible for pressures less than
  //    atmospheric.  To do that, we keep the accumulation derivative
  //    non-zero, calculating dWC_bar / dh_bar, where bar indicates (p -
  //    p_atm), not max(p - p_atm,0).  Note that this operator is in h
  //    coordinates, not p coordinates, as the diffusion operator is applied
  //    to h.
  //
  // -- update dh_bar / dp
  S_->GetEvaluator(pd_bar_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
  auto dh_dp = S_->GetDerivativePtr<CompositeVector>(pd_bar_key_, tag_next_, key_, tag_next_);

  // -- update the accumulation derivatives
  S_->GetEvaluator(wc_bar_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
  auto dwc_dp = S_->GetDerivativePtr<CompositeVector>(wc_bar_key_, tag_next_, key_, tag_next_);
  db_->WriteVector("    dwc_dp", dwc_dp.ptr());
  db_->WriteVector("    dh_dp", dh_dp.ptr());

  CompositeVector dwc_dh(dwc_dp->Map());
  dwc_dh.ReciprocalMultiply(1. / h, *dh_dp, *dwc_dp, 0.);
  preconditioner_acc_->AddAccumulationTerm(dwc_dh, "cell");

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(h);

  // 3. Assemble and precompute the Schur complement for inversion.
  // 3.a: Patch up BCs in the case of zero conductivity
  FixBCsForPrecon_(tag_next_);
  preconditioner_diff_->ApplyBCs(true, true, true);

  // 3.d: Rescale to use as a pressure matrix if used in a coupler
  //  if (coupled_to_subsurface_via_head_ || coupled_to_subsurface_via_flux_) {
  if (precon_scaled_) {
    // Scale Spp by dh/dp (h, NOT h_bar), clobbering rows with p < p_atm
    S_->GetEvaluator(pd_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
    Teuchos::RCP<const CompositeVector> dh0_dp =
      S_->GetDerivativePtr<CompositeVector>(pd_key_, tag_next_, key_, tag_next_);
    preconditioner_->Rescale(*dh0_dp);

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Right scaling TPFA" << std::endl;
    db_->WriteVector("    dh_dp", dh0_dp.ptr());
  }

  // increment the iterator count
  iter_++;
};

} // namespace Flow
} // namespace Amanzi
