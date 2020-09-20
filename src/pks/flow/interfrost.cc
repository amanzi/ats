/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */

#include "Op.hh"
#include "interfrost.hh"

namespace Amanzi {
namespace Flow {

// -- accumulation term
void
Interfrost::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  Permafrost::AddAccumulation_(g);
  double dt = S_next_->time() - S_inter_->time();

  // addition dp/dt part
  S_next_->GetFieldEvaluator("DThetaDp_coef")->HasFieldChanged(S_next_.ptr(), name_);
  S_next_->GetFieldEvaluator(key_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(key_)->HasFieldChanged(S_inter_.ptr(), name_);

  const Epetra_MultiVector& pres1 = *S_next_->GetFieldData(key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pres0 = *S_inter_->GetFieldData(key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv1 = *S_next_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);  
  const Epetra_MultiVector& dThdp = *S_next_->GetFieldData("DThetaDp_coef")
      ->ViewComponent("cell",false);

  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);
  for (int c=0; c!=g_c.MyLength(); ++c) {
    g_c[0][c] += cv1[0][c] * dThdp[0][c] * (pres1[0][c]-pres0[0][c]) / dt;
  }
}

void
Interfrost::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
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
  AMANZI_ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);
  //PKDefaultBase::solution_to_state(*up, S_next_);

  // update the rel perm according to the scheme of choice, also upwind derivatives of rel perm
  UpdatePermeabilityData_(S_next_.ptr());
  UpdatePermeabilityDerivativeData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData(uw_coef_key_);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    const Epetra_MultiVector& kr = *rel_perm->ViewComponent("face",false);
    double min_kr = 1.e6;
    int min_kr_lid = -1;
    for (int f=0; f!=kr.MyLength(); ++f) {
      if (kr[0][f] < min_kr) {
        min_kr = kr[0][f];
        min_kr_lid = f;
      }
    }
    AMANZI_ASSERT(min_kr_lid >= 0);

    ENorm_t global_min_kr;
    ENorm_t local_min_kr;
    local_min_kr.value = min_kr;
    local_min_kr.gid = kr.Map().GID(min_kr_lid);
    AMANZI_ASSERT(local_min_kr.gid >= 0);

    MPI_Allreduce(&local_min_kr, &global_min_kr, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    *vo_->os() << "Min Kr[face=" << global_min_kr.gid << "] = " << global_min_kr.value << std::endl;
  }

  Teuchos::RCP<CompositeVector> rel_perm_modified =
      Teuchos::rcp(new CompositeVector(*rel_perm));
  *rel_perm_modified = *rel_perm;

  {
    Epetra_MultiVector& rel_perm_mod_f = *rel_perm_modified->ViewComponent("face",false);
    unsigned int nfaces = rel_perm_mod_f.MyLength();
    for (unsigned int f=0; f!=nfaces; ++f) {
      rel_perm_mod_f[0][f] = std::max(rel_perm_mod_f[0][f], 1.e-18);
    }
  }      
  
  // Update the preconditioner with darcy and gravity fluxes
  preconditioner_->Init();

  S_next_->GetFieldEvaluator(mass_dens_key_)->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(mass_dens_key_);
  preconditioner_diff_->SetDensity(rho);

  Key dkrdp_key = Keys::getDerivKey(uw_coef_key_, key_);
  Teuchos::RCP<const CompositeVector> dkrdp = S_next_->GetFieldData(dkrdp_key);
   preconditioner_diff_->SetScalarCoefficient(rel_perm_modified, dkrdp);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  Teuchos::RCP<CompositeVector> flux = S_next_->GetFieldData(flux_key_, name_);
  preconditioner_diff_->UpdateFlux(up->Data().ptr(), flux.ptr());
  preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), Teuchos::null);
  
  // Update the preconditioner with accumulation terms.
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator(conserved_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);

  // -- get the accumulation deriv
  Key dwc_dp_key = Keys::getDerivKey(conserved_key_, key_);
  const Epetra_MultiVector& dwc_dp =
      *S_next_->GetFieldData(dwc_dp_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& pres =
      *S_next_->GetFieldData(key_)->ViewComponent("cell",false);

#if DEBUG_FLAG
  db_->WriteVector("    dwc_dp", S_next_->GetFieldData(dwc_dp_key).ptr());
#endif

  // -- and the extra interfrost deriv
  S_next_->GetFieldEvaluator("DThetaDp_coef")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& dThdp_coef =
      *S_next_->GetFieldData("DThetaDp_coef")->ViewComponent("cell",false);
  const Epetra_MultiVector& d_dThdp_coef_dp =
      *S_next_->GetFieldData("dDThetaDp_coef_dpressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pres0 =
      *S_inter_->GetFieldData(key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& cv =
      *S_next_->GetFieldData("cell_volume")->ViewComponent("cell",false);
  
  // -- update the cell-cell block

  auto& Acc_cells = *preconditioner_acc_->local_op(0)->diag;

  unsigned int ncells = dwc_dp.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    Acc_cells[0][c] += dwc_dp[0][c] / h + cv[0][c] * dThdp_coef[0][c] / h
        + d_dThdp_coef_dp[0][c] * cv[0][c] * (pres[0][c] - pres0[0][c])/h;
  }
  
  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(S_next_.ptr(), h);
  
  // -- apply BCs
  preconditioner_diff_->ApplyBCs(true, true, true);
}


// Create of physical evaluators.
void
Interfrost::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  Permafrost::SetupPhysicalEvaluators_(S);

  S->RequireField("DThetaDp_coef")
      ->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("DThetaDp_coef");
}



} // namespace Flow
} // namespace Amanzi
