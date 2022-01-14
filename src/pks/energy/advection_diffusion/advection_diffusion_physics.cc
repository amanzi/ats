/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Physics methods of the advection-diffusion equation.

Solves:

dT/dt  +  u dot grad(T) = div K grad T

------------------------------------------------------------------------- */

#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {


// dT/dt portion of the residual function
void AdvectionDiffusion::AddAccumulation_(Teuchos::RCP<CompositeVector> g) {
  S_next_->GetEvaluator("temperature")->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetEvaluator("temperature")->HasFieldChanged(S_inter_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> temp0 =
    S_inter_->GetPtr<CompositeVector>("temperature");
  Teuchos::RCP<const CompositeVector> temp1 =
    S_next_->GetPtr<CompositeVector>("temperature");

  Teuchos::RCP<const CompositeVector> cv0 =
    S_inter_->GetPtr<CompositeVector>("cell_volume");
  Teuchos::RCP<const CompositeVector> cv1 =
    S_next_->GetPtr<CompositeVector>("cell_volume");

  double dt = S_->get_time(tag_next_) - S_->get_time(tag_inter_);
  AMANZI_ASSERT(dt > 0.);

  //  --   g <-- g - (cv*h)_t0/dt
  g->ViewComponent("cell",false)->Multiply(-1./dt,
          *cv0->ViewComponent("cell",false), *temp0->ViewComponent("cell",false), 1.);
  //  --   g <-- g + (cv*h)_t1/dt
  g->ViewComponent("cell",false)->Multiply(1./dt,
          *cv1->ViewComponent("cell",false), *temp1->ViewComponent("cell",false), 1.);
};



// u dot grad T portion of the residual function
void AdvectionDiffusion::AddAdvection_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> g, bool negate) {

  // set up the operator
  Teuchos::RCP<const CompositeVector> mass_flux = S->GetPtr<CompositeVector>("mass_flux");
  matrix_adv_->global_operator()->Init();
  matrix_adv_->Setup(*mass_flux);
  matrix_adv_->SetBCs(bc_, bc_);
  matrix_adv_->UpdateMatrices(mass_flux.ptr());
  matrix_adv_->ApplyBCs(false, true, false);


  // apply
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>("temperature");
  matrix_adv_->global_operator()->ComputeNegativeResidual(*temp, *g, false);
};


// - div K grad T part of the residual
void AdvectionDiffusion::ApplyDiffusion_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> g) {
  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S->GetPtr<CompositeVector>("thermal_conductivity");

  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(thermal_conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // calculate the div-grad operator, apply it to temperature, and add to residual
  Teuchos::RCP<const CompositeVector> temp =
    S->GetPtr<CompositeVector>("temperature");


  // finish assembly of the stiffness matrix
  matrix_diff_->ApplyBCs(true, true, true);

  // calculate the residual
  matrix_diff_->global_operator()->ComputeNegativeResidual(*temp, *g);
};

} //namespace Energy
} //namespace Amanzi
