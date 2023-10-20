/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Op.hh"
#include "interfrost.hh"

namespace Amanzi {
namespace Flow {

// -- accumulation term
void
Interfrost::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  Permafrost::AddAccumulation_(g);
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  // addition dp/dt part
  S_->GetEvaluator("DThetaDp_coef", tag_next_).Update(*S_, name_);
  S_->GetEvaluator(key_, tag_next_).Update(*S_, name_);
  S_->GetEvaluator(key_, tag_current_).Update(*S_, name_);

  const Epetra_MultiVector& pres1 =
    *S_->Get<CompositeVector>(key_, tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& pres0 =
    *S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);
  const Epetra_MultiVector& cv1 =
    *S_->Get<CompositeVector>("cell_volume", tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& dThdp =
    *S_->Get<CompositeVector>("DThetaDp_coef", tag_next_).ViewComponent("cell", false);

  Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);
  for (int c = 0; c != g_c.MyLength(); ++c) {
    g_c[0][c] += cv1[0][c] * dThdp[0][c] * (pres1[0][c] - pres0[0][c]) / dt;
  }
}

void
Interfrost::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon update at t = " << t << std::endl;

  // Recreate mass matrices
  if (!deform_key_.empty() &&
      S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " precon"))
    preconditioner_diff_->SetTensorCoefficient(K_);

  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);

  // update the rel perm according to the scheme of choice, also upwind derivatives of rel perm
  UpdatePermeabilityData_(tag_next_);
  UpdatePermeabilityDerivativeData_(tag_next_);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);

  // modify the rel perm according to the interfrost model -- HACK!
  // NOTE: this may not be necessary anymore, but is left here for legacy reasons
  const CompositeVector& rel_perm = S_->Get<CompositeVector>(uw_coef_key_, tag_next_);
  auto rel_perm_modified = Teuchos::rcp(new CompositeVector(rel_perm));
  *rel_perm_modified = rel_perm;
  {
    Epetra_MultiVector& rel_perm_mod_f = *rel_perm_modified->ViewComponent("face", false);
    unsigned int nfaces = rel_perm_mod_f.MyLength();
    for (unsigned int f = 0; f != nfaces; ++f) {
      rel_perm_mod_f[0][f] = std::max(rel_perm_mod_f[0][f], 1.e-18);
    }
  }

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
  preconditioner_diff_->SetScalarCoefficient(rel_perm_modified, dkrdp);

  // -- fill local matrices
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  preconditioner_diff_->ApplyBCs(true, true, true);

  // -- update with Jacobian terms
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

  // -- and the extra interfrost deriv
  S_->GetEvaluator("DThetaDp_coef", tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
  const Epetra_MultiVector& dwc_dp_vec = *dwc_dp->ViewComponent("cell", false);
  const Epetra_MultiVector& dThdp_coef =
    *S_->Get<CompositeVector>("DThetaDp_coef", tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& d_dThdp_coef_dp =
    *S_->GetDerivative<CompositeVector>("dDThetaDp_coef", tag_next_, key_, tag_next_)
       .ViewComponent("cell", false);
  const Epetra_MultiVector& pres0 =
    *S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);
  const Epetra_MultiVector& pres1 =
    *S_->Get<CompositeVector>(key_, tag_next_).ViewComponent("cell", false);
  const Epetra_MultiVector& cv =
    *S_->Get<CompositeVector>("cell_volume", tag_next_).ViewComponent("cell", false);

  // -- update the cell-cell block
  auto& Acc_cells = *preconditioner_acc_->local_op(0)->diag;
  unsigned int ncells = dwc_dp_vec.MyLength();
  for (unsigned int c = 0; c != ncells; ++c) {
    Acc_cells[0][c] += dwc_dp_vec[0][c] / h + cv[0][c] * dThdp_coef[0][c] / h +
                       d_dThdp_coef_dp[0][c] * cv[0][c] * (pres1[0][c] - pres0[0][c]) / h;
  }

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(h);

  // increment the iterator count
  iter_++;
}


// Create of physical evaluators.
void
Interfrost::SetupPhysicalEvaluators_()
{
  Permafrost::SetupPhysicalEvaluators_();

  // in addition, require DThetaDP_coef, the specific storage term in Interfrost model
  S_->Require<CompositeVector, CompositeVectorSpace>("DThetaDp_coef", tag_next_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator("DThetaDp_coef", tag_next_);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    "DThetaDp_coef", tag_next_, key_, tag_next_);
}


} // namespace Flow
} // namespace Amanzi
