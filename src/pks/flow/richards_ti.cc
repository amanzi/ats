/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "Op.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Richards::FunctionalResidual(double t_old,
                   double t_new,
                   Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new,
                   Teuchos::RCP<TreeVector> g)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  double h = t_new - t_old;

  //--  AMANZI_ASSERT(std::abs(S_->get_time(tag_inter_) - t_old) < 1.e-4*h);
  //-- AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4*h);

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, tag_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

  if (dynamic_mesh_) matrix_diff_->SetTensorCoefficient(K_);

  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << h << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("p_old"); vnames.push_back("p_new");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_->GetPtr<CompositeVector>(key_, tag_current_).ptr()); vecs.push_back(u.ptr());
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
  // if (vapor_diffusion_) AddVaporDiffusionResidual_(tag_next_, res.ptr());

  // dump s_old, s_new
  vnames[0] = "sl_old"; vnames[1] = "sl_new";
  vecs[0] = S_->GetPtr<CompositeVector>(sat_key_, tag_current_).ptr();
  vecs[1] = S_->GetPtr<CompositeVector>(sat_key_, tag_next_).ptr();

  if (S_->HasRecordSet(sat_ice_key_)) {
    vnames.push_back("si_old");
    vnames.push_back("si_new");
    vecs.push_back(S_->GetPtr<CompositeVector>(Keys::getKey(domain_,"saturation_ice"), tag_current_).ptr());
    vecs.push_back(S_->GetPtr<CompositeVector>(Keys::getKey(domain_,"saturation_ice"), tag_next_).ptr());
  }
  vnames.push_back("poro");
  vecs.push_back(S_->GetPtr<CompositeVector>(Keys::getKey(domain_,"porosity"), tag_next_).ptr());
  vnames.push_back("perm_K");
  vecs.push_back(S_->GetPtr<CompositeVector>(Keys::getKey(domain_,"permeability"), tag_next_).ptr());
  vnames.push_back("k_rel");
  vecs.push_back(S_->GetPtr<CompositeVector>(coef_key_, tag_next_).ptr());
  vnames.push_back("wind");
  vecs.push_back(S_->GetPtr<CompositeVector>(flux_dir_key_, tag_next_).ptr());
  vnames.push_back("uw_k_rel");
  vecs.push_back(S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_).ptr());
  vnames.push_back("flux");
  vecs.push_back(S_->GetPtr<CompositeVector>(flux_key_, tag_next_).ptr());
  db_->WriteVectors(vnames,vecs,true);

  db_->WriteVector("res (diff)", res.ptr(), true);

  // accumulation term
  AddAccumulation_(res.ptr());
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
int Richards::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

  db_->WriteVector("p_res", u->Data().ptr(), true);

  // Apply the preconditioner
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

  db_->WriteVector("PC*p_res", Pu->Data().ptr(), true);
  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Richards::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // Recreate mass matrices
  if (dynamic_mesh_) {
    matrix_diff_->SetTensorCoefficient(K_);
    preconditioner_diff_->SetTensorCoefficient(K_);
  }

  // update state with the solution up.
  if (std::abs(t - iter_counter_time_)/t > 1.e-4) {
    iter_ = 0;
    iter_counter_time_ = t;
  }
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4*t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);

  // update the rel perm according to the scheme of choice, also upwind derivatives of rel perm
  UpdatePermeabilityData_(tag_next_);
  if (jacobian_ && iter_ >= jacobian_lag_) UpdatePermeabilityDerivativeData_(tag_next_);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);

  Teuchos::RCP<const CompositeVector> rel_perm =
    S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_);

  // Update the preconditioner with darcy and gravity fluxes
  preconditioner_->Init();

  // gravity fluxes
  S_->GetEvaluator(mass_dens_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag_next_);
  preconditioner_diff_->SetDensity(rho);

  // jacobian term
  Teuchos::RCP<const CompositeVector> dkrdp = Teuchos::null;
  if (jacobian_ && iter_ >= jacobian_lag_) {
    if (!duw_coef_key_.empty()) {
      dkrdp = S_->GetPtr<CompositeVector>(duw_coef_key_, tag_next_);
    } else {
      dkrdp = S_->GetPtr<CompositeVector>(dcoef_key_, tag_next_);
    }
  }

  // create local matrices
  preconditioner_diff_->SetScalarCoefficient(rel_perm, dkrdp);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  preconditioner_diff_->ApplyBCs(true, true, true);

  if (jacobian_ && iter_ >= jacobian_lag_) {// && preconditioner_->RangeMap().HasComponent("face")) {
    Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(flux_key_, tag_next_, name_);
    preconditioner_diff_->UpdateFlux(up->Data().ptr(), flux.ptr());
    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), up->Data().ptr());
  }

  // Update the preconditioner with accumulation terms.
  // -- update the accumulation derivatives
  S_->GetEvaluator(conserved_key_, tag_next_).UpdateDerivative(*S_, name_,
          key_, tag_next_);

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> dwc_dp =
    S_->GetDerivativePtr<CompositeVector>(conserved_key_, tag_next_,
            key_, tag_next_);
  db_->WriteVector("    dwc_dp", dwc_dp.ptr());

  // -- update the cell-cell block  CompositeVector du(S_->Get<CompositeVector>(dwc_dp_key, tag_next_).Map());
  preconditioner_acc_->AddAccumulationTerm(*dwc_dp, h, "cell", false);

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(h);

  // increment the iterator count
  iter_++;
};


}  // namespace Flow
}  // namespace Amanzi



