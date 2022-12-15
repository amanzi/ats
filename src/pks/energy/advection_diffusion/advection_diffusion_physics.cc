/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Physics methods of the advection-diffusion equation.

Solves:

dT/dt  +  u dot grad(T) = div K grad T

------------------------------------------------------------------------- */

#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {


// dT/dt portion of the residual function
void
AdvectionDiffusion::AddAccumulation_(Teuchos::RCP<CompositeVector> g)
{
  S_next_->GetFieldEvaluator("temperature")->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator("temperature")->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> temp0 = S_inter_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> temp1 = S_next_->GetFieldData("temperature");

  Teuchos::RCP<const CompositeVector> cv0 = S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cv1 = S_next_->GetFieldData("cell_volume");

  double dt = S_next_->time() - S_inter_->time();
  AMANZI_ASSERT(dt > 0.);

  //  --   g <-- g - (cv*h)_t0/dt
  g->ViewComponent("cell", false)
    ->Multiply(
      -1. / dt, *cv0->ViewComponent("cell", false), *temp0->ViewComponent("cell", false), 1.);
  //  --   g <-- g + (cv*h)_t1/dt
  g->ViewComponent("cell", false)
    ->Multiply(
      1. / dt, *cv1->ViewComponent("cell", false), *temp1->ViewComponent("cell", false), 1.);
};


// u dot grad T portion of the residual function
void
AdvectionDiffusion::AddAdvection_(const Teuchos::RCP<State> S,
                                  const Teuchos::RCP<CompositeVector> g,
                                  bool negate)
{
  // set up the operator
  Teuchos::RCP<const CompositeVector> water_flux = S->GetFieldData("water_flux");
  matrix_adv_->global_operator()->Init();
  matrix_adv_->Setup(*water_flux);
  matrix_adv_->SetBCs(bc_, bc_);
  matrix_adv_->UpdateMatrices(water_flux.ptr());
  matrix_adv_->ApplyBCs(false, true, false);


  // apply
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  matrix_adv_->global_operator()->ComputeNegativeResidual(*temp, *g, false);
};


// - div K grad T part of the residual
void
AdvectionDiffusion::ApplyDiffusion_(const Teuchos::RCP<State> S,
                                    const Teuchos::RCP<CompositeVector> g)
{
  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity");

  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(thermal_conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // calculate the div-grad operator, apply it to temperature, and add to residual
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");


  // finish assembly of the stiffness matrix
  matrix_diff_->ApplyBCs(true, true, true);

  // calculate the residual
  matrix_diff_->global_operator()->ComputeNegativeResidual(*temp, *g);
};

} //namespace Energy
} //namespace Amanzi
