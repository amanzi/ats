/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "overland_pressure.hh"
#include "Op.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_RES_FLAG 0

// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void OverlandPressureFlow::FunctionalResidual( double t_old,
                        double t_new,
                        Teuchos::RCP<TreeVector> u_old,
                        Teuchos::RCP<TreeVector> u_new,
                        Teuchos::RCP<TreeVector> g )
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  niter_++;

  // bookkeeping
  double h = t_new - t_old;

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, S_next_);

  Teuchos::RCP<CompositeVector> u = u_new->Data();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << h << std::endl;

    // double minv, maxv;
    // u->MinValue(&minv);
    // u->MaxValue(&maxv);
    // *vo_->os() << "  min,max_u = " << minv << ", " << maxv << std::endl;
  }

  // unnecessary here if not debeugging, but doesn't hurt either
  S_next_->GetFieldEvaluator(potential_key_)->HasFieldChanged(S_next_.ptr(), name_);

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("p_old");
  vnames.push_back("p_new");
  vnames.push_back("z");
  vnames.push_back("h_old");
  vnames.push_back("h_new");
  vnames.push_back("h+z");
  if (plist_->isSublist("overland conductivity subgrid evaluator")) {
    vnames.push_back("pd - dd");
    vnames.push_back("frac_cond");
  }

  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData(key_).ptr());
  vecs.push_back(u.ptr());

  vecs.push_back(S_inter_->GetFieldData(elev_key_).ptr());
  vecs.push_back(S_inter_->GetFieldData(pd_key_).ptr());
  vecs.push_back(S_next_->GetFieldData(pd_key_).ptr());
  vecs.push_back(S_next_->GetFieldData(potential_key_).ptr());

  if (plist_->isSublist("overland conductivity subgrid evaluator")) {
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"mobile_depth")).ptr());
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"fractional_conductance")).ptr());
  }
  db_->WriteVectors(vnames, vecs, true);

  // update boundary conditions
  bc_head_->Compute(S_next_->time());
  bc_pressure_->Compute(S_next_->time());
  bc_zero_gradient_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  bc_level_->Compute(S_next_->time());
  bc_seepage_head_->Compute(S_next_->time());
  bc_seepage_pressure_->Compute(S_next_->time());
  bc_critical_depth_->Compute(S_next_->time());
  bc_dynamic_->Compute(S_next_->time());
  bc_tidal_->Compute(S_next_->time());

  UpdateBoundaryConditions_(S_next_.ptr());

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

  db_->WriteBoundaryConditions(bc_markers(), bc_values());
  if (S_next_->HasField(Keys::getKey(domain_,"unfrozen_fraction"))) {
    Key uf_key = Keys::getKey(domain_,"unfrozen_fraction");
    S_next_->GetFieldEvaluator(uf_key)->HasFieldChanged(S_next_.ptr(), name_);
    vnames.resize(2);
    vecs.resize(2);
    vnames[0] = "uf_frac_old";
    vnames[1] = "uf_frac_new";
    vecs[0] = S_inter_->GetFieldData(uf_key).ptr();
    vecs[1] = S_next_->GetFieldData(uf_key).ptr();
    db_->WriteVectors(vnames, vecs, true);
  }
  db_->WriteVector("uw_dir", S_next_->GetFieldData(flux_dir_key_).ptr(), true);
  db_->WriteVector("k_s", S_next_->GetFieldData(cond_key_).ptr(), true);
  db_->WriteVector("k_s_uw", S_next_->GetFieldData(uw_cond_key_).ptr(), true);
  db_->WriteVector("q_s", S_next_->GetFieldData(flux_key_).ptr(), true);
  db_->WriteVector("res (diff)", res.ptr(), true);

  // accumulation term
  AddAccumulation_(res.ptr());
  db_->WriteVector("res (acc)", res.ptr(), true);

  // add rhs load value
  AddSourceTerms_(res.ptr());
  db_->WriteVector("res (src)", res.ptr(), true);

#if DEBUG_RES_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    double minv, maxv;
    res->MinValue(&minv);
    res->MaxValue(&maxv);
    *vo_->os() << "  min,max_res = " << minv << ", " << maxv << std::endl;
  }

  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << "surface-flow_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::cout << "WRITING FILE: " << namestream.str() << std::endl;
    Teuchos::RCP<const CompositeVector> depth= S_next_->GetFieldData(pd_key_);

    std::stringstream solnstream;
    solnstream << "surface-flow_solution_" << niter_;
    std::cout << "WRITING FILE: " << solnstream.str() << std::endl;
    *S_next_->GetFieldData(solnstream.str(),name_) = *depth;
  }
#endif

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int OverlandPressureFlow::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;
  AMANZI_ASSERT(!precon_scaled_); // otherwise this factor was built into the matrix

  // apply the preconditioner
  db_->WriteVector("h_res", u->Data().ptr(), true);
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
  db_->WriteVector("PC*h_res (h-coords)", Pu->Data().ptr(), true);

  // tack on the variable change
  const Epetra_MultiVector& dh_dp =
    *S_next_->GetFieldData(Keys::getDerivKey(pd_bar_key_,key_))->ViewComponent("cell",false);
  Epetra_MultiVector& Pu_c = *Pu->Data()->ViewComponent("cell",false);

  unsigned int ncells = Pu_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    Pu_c[0][c] /= dh_dp[0][c];
  }

  db_->WriteVector("PC*h_res (p-coords)", Pu->Data().ptr(), true);
  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandPressureFlow::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  if (std::abs(t - iter_counter_time_)/t > 1.e-4) {
    iter_ = 0;
    iter_counter_time_ = t;
  }

  AMANZI_ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);

  // calculating the operator is done in 3 steps:
  // 1. Diffusion components

  // 1.a: Pre-assembly updates.
  // -- update boundary condition markers, which set the BC type
  UpdateBoundaryConditions_(S_next_.ptr());

  // -- update the rel perm according to the boundary info and upwinding
  // -- scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());
  if (jacobian_ && iter_ >= jacobian_lag_) UpdatePermeabilityDerivativeData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData(uw_cond_key_);

  Teuchos::RCP<const CompositeVector> dcond = Teuchos::null;
  if (jacobian_ && iter_ >= jacobian_lag_) {
    if (preconditioner_->RangeMap().HasComponent("face")) {
      dcond = S_next_->GetFieldData(duw_cond_key_);
    } else {
      dcond = S_next_->GetFieldData(dcond_key_);
    }
  }

  // 1.b: Create all local matrices.
  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(cond, dcond);

  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  if (jacobian_ && iter_ >= jacobian_lag_) {
    Teuchos::RCP<const CompositeVector> pres_elev = Teuchos::null;
    Teuchos::RCP<CompositeVector> flux = Teuchos::null;
    if (preconditioner_->RangeMap().HasComponent("face")) {
      flux = S_next_->GetFieldData(flux_key_, name_);
      preconditioner_diff_->UpdateFlux(pres_elev.ptr(), flux.ptr());
    } else {
      S_next_->GetFieldEvaluator(potential_key_)->HasFieldChanged(S_next_.ptr(), name_);
      pres_elev = S_next_->GetFieldData(potential_key_);
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
  S_next_->GetFieldEvaluator(pd_bar_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  auto dh_dp = S_next_->GetFieldData(Keys::getDerivKey(pd_bar_key_, key_));

  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator(wc_bar_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  auto dwc_dp = S_next_->GetFieldData(Keys::getDerivKey(wc_bar_key_, key_));
  db_->WriteVector("    dwc_dp", dwc_dp.ptr());
  db_->WriteVector("    dh_dp", dh_dp.ptr());

  CompositeVector dwc_dh(dwc_dp->Map());
  dwc_dh.ReciprocalMultiply(1./h, *dh_dp, *dwc_dp, 0.);
  preconditioner_acc_->AddAccumulationTerm(dwc_dh, "cell");

  // Why is this turned off? #60 --etc
  // // -- update the source term derivatives
  // if (S_next_->GetFieldEvaluator(source_key_)->IsDependency(S_next_.ptr(), key_)) {
  //   S_next_->GetFieldEvaluator(source_key_)
  //       ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  //   Key dkey = Keys::getDerivKey(source_key_,key_);
  //   const Epetra_MultiVector& dq_dp = *S_next_->GetFieldData(dkey)
  //       ->ViewComponent("cell",false);

  //   const Epetra_MultiVector& cv =
  //       *S_next_->GetFieldData("surface-cell_volume")->ViewComponent("cell",false);

  //   if (source_in_meters_) {
  //     // External source term is in [m water / s], not in [mols / s], so a
  //     // density is required.  This density should be upwinded.
  //     S_next_->GetFieldEvaluator("surface-molar_density_liquid")
  //         ->HasFieldChanged(S_next_.ptr(), name_);
  //     S_next_->GetFieldEvaluator("surface-source_molar_density")
  //         ->HasFieldChanged(S_next_.ptr(), name_);
  //     const Epetra_MultiVector& nliq1 =
  //         *S_next_->GetFieldData("surface-molar_density_liquid")
  //         ->ViewComponent("cell",false);
  //     const Epetra_MultiVector& nliq1_s =
  //       *S_next_->GetFieldData("surface-source_molar_density")
  //         ->ViewComponent("cell",false);
  //     const Epetra_MultiVector& q = *S_next_->GetFieldData(source_key_)
  //         ->ViewComponent("cell",false);

  //     for (int c=0; c!=cv.MyLength(); ++c) {
  //       double s1 = q[0][c] > 0. ? dq_dp[0][c] * nliq1_s[0][c] : dq_dp[0][c] * nliq1[0][c];
  //       Acc_cells[c] -= cv[0][c] * s1 / dh_dp[0][c];
  //     }
  //   } else {
  //     for (int c=0; c!=cv.MyLength(); ++c) {
  //       Acc_cells[c] -= cv[0][c] * dq_dp[0][c] / dh_dp[0][c];
  //     }
  //   }
  // }


  // 3. Assemble and precompute the Schur complement for inversion.
  // 3.a: Patch up BCs in the case of zero conductivity
  FixBCsForPrecon_(S_next_.ptr());

  preconditioner_diff_->ApplyBCs(true, true, true);

  // 3.d: Rescale to use as a pressure matrix if used in a coupler
  //  if (coupled_to_subsurface_via_head_ || coupled_to_subsurface_via_flux_) {
  if (precon_scaled_) {
    // Scale Spp by dh/dp (h, NOT h_bar), clobbering rows with p < p_atm
    std::string pd_deriv_key = Keys::getDerivKey(pd_key_,key_);

    S_next_->GetFieldEvaluator(pd_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    Teuchos::RCP<const CompositeVector> dh0_dp = S_next_->GetFieldData(pd_deriv_key);
    const Epetra_MultiVector& dh0_dp_c = *dh0_dp->ViewComponent("cell",false);

    preconditioner_->Rescale(*dh0_dp);

    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  Right scaling TPFA" << std::endl;
    db_->WriteVector("    dh_dp", dh0_dp.ptr());
  }

  // increment the iterator count
  iter_++;
};

}  // namespace Flow
}  // namespace Amanzi
