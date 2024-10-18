/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "Epetra_Vector.h"
#include "advection_diffusion.hh"
#include "Op.hh"
#include "EpetraExt_RowMatrixOut.h"

namespace Amanzi {
namespace Energy {

// AdvectionDiffusion is a BDFFnBase
// computes the non-linear functional g = g(t,u,udot)
void
AdvectionDiffusion::FunctionalResidual(double t_old,
                                       double t_new,
                                       Teuchos::RCP<const TreeVector> u_old,
                                       Teuchos::RCP<TreeVector> u_new,
                                       Teuchos::RCP<TreeVector> g)
{
  // pointer-copy temperature into states and update any auxilary data
  Solution_to_State(*u_new, S_next_);

  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  Teuchos::RCP<CompositeVector> u = u_new->Data();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Residual calculation:" << std::endl;
    *vo_->os() << "  u: " << (*u)("cell", 0) << std::endl;
  }

  // get access to the solution
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_, res);
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  res (after diffusion): " << (*res)("cell", 0) << "," << (*res)("cell", 19)
               << std::endl;

  // accumulation term
  AddAccumulation_(res);
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  res (after accumulation): " << (*res)("cell", 0) << "," << (*res)("cell", 19)
               << std::endl;

  // advection term, explicit
  if (implicit_advection_) {
    AddAdvection_(S_next_, res, true);
  } else {
    AddAdvection_(S_inter_, res, true);
  }
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  res (after advection): " << (*res)("cell", 0) << "," << (*res)("cell", 19)
               << std::endl;
};

// applies preconditioner to u and returns the result in Pu
int
AdvectionDiffusion::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                        Teuchos::RCP<TreeVector> Pu)
{
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Precon application:" << std::endl;
    *vo_->os() << "  u: " << (*u->Data())("cell", 0);
    if (u->Data()->HasComponent("face")) *vo_->os() << "  f: " << (*u->Data())("face", 80);
    *vo_->os() << std::endl;
  }

  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "  Pu: " << (*Pu->Data())("cell", 0);
    if (Pu->Data()->HasComponent("face")) *vo_->os() << "  f: " << (*Pu->Data())("face", 80);
    *vo_->os() << std::endl;
  }

  return (ierr > 0) ? 0 : 1;
};


// updates the preconditioner
void
AdvectionDiffusion::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  AMANZI_ASSERT(std::abs(S_next_->time() - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // div K_e grad u
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S_next_->GetFieldData("thermal_conductivity");
  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(thermal_conductivity, Teuchos::null);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // update with accumulation terms
  Teuchos::RCP<const CompositeVector> cell_volume = S_next_->GetFieldData("cell_volume");

  CompositeVector du(cell_volume->Map());
  du.PutScalar(1.);
  preconditioner_acc_->AddAccumulationTerm(du, h, "cell");

  // update with advection terms
  if (implicit_advection_) {
    Teuchos::RCP<const CompositeVector> water_flux = S_next_->GetFieldData("water_flux");
    preconditioner_adv_->Setup(*water_flux);
    preconditioner_adv_->SetBCs(bc_, bc_);
    preconditioner_adv_->UpdateMatrices(water_flux.ptr());
    preconditioner_adv_->ApplyBCs(false, true, false);
  }

  // assemble and create PC
  preconditioner_diff_->ApplyBCs(true, true, true);
  preconditioner_->AssembleMatrix();
  preconditioner_->UpdatePreconditioner();
};


} // namespace Energy
} // namespace Amanzi
