/* -*-  mode++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Solves:

de/dt + q dot grad h = div Ke grad T + S?
------------------------------------------------------------------------- */

#include "BoundaryFunction.hh"
#include "Evaluator.hh"
#include "Op.hh"
#include "energy_interfrost.hh"

namespace Amanzi {
namespace Energy {


void
InterfrostEnergy::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  ThreePhase::SetupPhysicalEvaluators_(S);

  S->Require<CompositeVector,CompositeVectorSpace>("DEnergyDT_coef", Tags::NEXT)
      ->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator("DEnergyDT_coef");

}


// -------------------------------------------------------------
// Accumulation of energy term de/dt
// -------------------------------------------------------------
void
InterfrostEnergy::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_inter_);

  // update the energy at both the old and new times.
  S_next_->GetEvaluator("DEnergyDT_coef")->HasFieldChanged(S_next_.ptr(), name_);
  S_next_->GetEvaluator(key_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetEvaluator(key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // get the energy at each time
  const Epetra_MultiVector& cv = *S_next_->Get<CompositeVector>(cell_vol_key_).ViewComponent("cell",false);
  const Epetra_MultiVector& dEdT_coef = *S_next_->Get<CompositeVector>("DEnergyDT_coef").ViewComponent("cell",false);
  const Epetra_MultiVector& T1 = *S_next_->Get<CompositeVector>(key_).ViewComponent("cell",false);
  const Epetra_MultiVector& T0 = *S_inter_->Get<CompositeVector>(key_).ViewComponent("cell",false);

  Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  for (int c=0; c!=g_c.MyLength(); ++c) {
    g_c[0][c] += cv[0][c] * dEdT_coef[0][c] * (T1[0][c] - T0[0][c]) / dt;
  }
};


void
InterfrostEnergy::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4*t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);
  Teuchos::RCP<const CompositeVector> temp = S_next_->GetPtr<CompositeVector>(key_);

  // update boundary conditions
  bc_temperature_->Compute(S_->get_time(tag_next_));
  bc_diff_flux_->Compute(S_->get_time(tag_next_));
  bc_flux_->Compute(S_->get_time(tag_next_));
  UpdateBoundaryConditions_(S_next_.ptr());

  // div K_e grad u
  UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetPtr<CompositeVector>(conductivity_key_);

  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, temp.ptr());

  // update with accumulation terms
  // -- update the accumulation derivatives, de/dT
  S_next_->GetEvaluator("DEnergyDT_coef")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& dcoef_dT = *S_next_->GetPtr<CompositeVector>("dDEnergyDT_coef_dtemperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& coef = *S_next_->GetPtr<CompositeVector>("DEnergyDT_coef")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S_next_->Get<CompositeVector>(cell_vol_key_).ViewComponent("cell",false);
  const Epetra_MultiVector& T1 = *S_next_->Get<CompositeVector>(key_).ViewComponent("cell",false);
  const Epetra_MultiVector& T0 = *S_inter_->Get<CompositeVector>(key_).ViewComponent("cell",false);


#if DEBUG_FLAG
  db_->WriteVector("    de_dT", S_next_->GetPtr<CompositeVector>(de_dT_key_).ptr());
#endif

  // -- get the matrices/rhs that need updating
  auto& Acc_cells = *preconditioner_acc_->local_op(0)->diag;

  // -- update the diagonal
  unsigned int ncells = T0.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    Acc_cells[0][c] += coef[0][c]*cv[0][c]/h + dcoef_dT[0][c] * cv[0][c] * (T1[0][c]-T0[0][c])/h;
  }

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(S_next_.ptr(), h);

  // update with advection terms
  if (implicit_advection_ && implicit_advection_in_pc_) {
    Teuchos::RCP<const CompositeVector> mass_flux = S_next_->GetPtr<CompositeVector>("mass_flux");
    S_next_->GetEvaluator(enthalpy_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    Teuchos::RCP<const CompositeVector> dhdT = S_next_->GetPtrW<CompositeVector>(Keys::getDerivKey(enthalpy_key_,key_));
    preconditioner_adv_->Setup(*mass_flux);
    preconditioner_adv_->UpdateMatrices(mass_flux.ptr(), dhdT.ptr());
    ApplyDirichletBCsToEnthalpy_(S_next_.ptr());
    preconditioner_adv_->ApplyBCs(false, true, false);
  }

  // Apply boundary conditions.
  preconditioner_diff_->ApplyBCs(true, true, true);
};



} // namespace Energy
} // namespace Amanzi
