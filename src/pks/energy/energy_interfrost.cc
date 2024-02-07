/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Solves:

de/dt + q dot grad h = div Ke grad T + S?
------------------------------------------------------------------------- */

#include "BoundaryFunction.hh"
#include "EvaluatorPrimary.hh"
#include "Op.hh"
#include "energy_interfrost.hh"

namespace Amanzi {
namespace Energy {


void
InterfrostEnergy::SetupPhysicalEvaluators_()
{
  ThreePhase::SetupPhysicalEvaluators_();

  S_->Require<CompositeVector, CompositeVectorSpace>("DEnergyDT_coef", tag_next_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator("DEnergyDT_coef", tag_next_);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    "DEnergyDT_coef", tag_next_, key_, tag_next_);
}


// -------------------------------------------------------------
// Accumulation of energy term de/dt
// -------------------------------------------------------------
void
InterfrostEnergy::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  // update the energy at both the old and new times.
  S_->GetEvaluator("DEnergyDT_coef", tag_next_).Update(*S_, name_);
  S_->GetEvaluator(key_, tag_next_).Update(*S_, name_);
  S_->GetEvaluator(key_, tag_current_).Update(*S_, name_);

  // get the energy at each time
  const auto& cv = *S_->Get<CompositeVector>(cell_vol_key_, tag_next_).ViewComponent("cell", false);
  const auto& dEdT_coef =
    *S_->Get<CompositeVector>("DEnergyDT_coef", tag_next_).ViewComponent("cell", false);
  const auto& T1 = *S_->Get<CompositeVector>(key_, tag_next_).ViewComponent("cell", false);
  const auto& T0 = *S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);

  auto& g_c = *g->ViewComponent("cell", false);

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  for (int c = 0; c != g_c.MyLength(); ++c) {
    g_c[0][c] += cv[0][c] * dEdT_coef[0][c] * (T1[0][c] - T0[0][c]) / dt;
  }
};


void
InterfrostEnergy::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);
  Teuchos::RCP<const CompositeVector> temp = S_->GetPtr<CompositeVector>(key_, tag_next_);

  // update thermal conductivity according to the boundary info and upwindin
  UpdateConductivityData_(tag_next_);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);

  // div K_e grad u
  Teuchos::RCP<const CompositeVector> conductivity =
    S_->GetPtr<CompositeVector>(conductivity_key_, tag_next_);

  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, temp.ptr());

  // update with accumulation terms
  // -- update the accumulation derivatives, de/dT
  S_->GetEvaluator("DEnergyDT_coef", tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
  const auto& dcoef_dT =
    *S_->GetDerivative<CompositeVector>("DEnergyDT_coef", tag_next_, key_, tag_next_)
       .ViewComponent("cell", false);
  const auto& coef =
    *S_->Get<CompositeVector>("DEnergyDT_coef", tag_next_).ViewComponent("cell", false);
  const auto& cv = *S_->Get<CompositeVector>(cell_vol_key_, tag_next_).ViewComponent("cell", false);
  const auto& T1 = *S_->Get<CompositeVector>(key_, tag_next_).ViewComponent("cell", false);
  const auto& T0 = *S_->Get<CompositeVector>(key_, tag_current_).ViewComponent("cell", false);

#if DEBUG_FLAG
  db_->WriteVector("    de_dT", S_->Get<CompositeVector>("DEnergyDT_coef", tag_next_).ptr());
#endif

  // -- get the matrices/rhs that need updating
  auto& Acc_cells = *preconditioner_acc_->local_op(0)->diag;

  // -- update the diagonal
  unsigned int ncells = T0.MyLength();
  for (unsigned int c = 0; c != ncells; ++c) {
    Acc_cells[0][c] +=
      coef[0][c] * cv[0][c] / h + dcoef_dT[0][c] * cv[0][c] * (T1[0][c] - T0[0][c]) / h;
  }

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(h);

  // update with advection terms
  if (implicit_advection_ && implicit_advection_in_pc_) {
    Teuchos::RCP<const CompositeVector> water_flux =
      S_->GetPtr<CompositeVector>("water_flux", tag_next_);
    S_->GetEvaluator(enthalpy_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
    const auto dhdT = S_->GetDerivativePtr<CompositeVector>(
      Keys::getDerivKey(enthalpy_key_, key_), tag_next_, key_, tag_next_);
    preconditioner_adv_->Setup(*water_flux);
    preconditioner_adv_->UpdateMatrices(water_flux.ptr(), dhdT.ptr());
    ApplyDirichletBCsToEnthalpy_(tag_next_);
    preconditioner_adv_->ApplyBCs(false, true, false);
  }

  // Apply boundary conditions.
  preconditioner_diff_->ApplyBCs(true, true, true);
};


} // namespace Energy
} // namespace Amanzi
