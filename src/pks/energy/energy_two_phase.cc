/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */


#include "eos_evaluator.hh"
#include "iem_evaluator.hh"
#include "thermal_conductivity_twophase_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "energy_bc_factory.hh"

#include "energy_two_phase.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------

TwoPhase::TwoPhase(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    EnergyBase(FElist, plist, S, solution) {}


// -------------------------------------------------------------
// Create the physical evaluators for energy and thermal conductivity
// -------------------------------------------------------------
void TwoPhase::SetupPhysicalEvaluators_() {
  // Get data and evaluators needed by the PK
  // -- energy, energy evaluator, and energy derivative
  S_->Require<CompositeVector,CompositeVectorSpace>(conserved_key_, tag_next_).SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(conserved_key_, tag_next_);
  S_->RequireDerivative<CompositeVector,CompositeVectorSpace>(conserved_key_, tag_next_, key_, tag_next_);

  // energy at the current time, where it is a copy evaluator
  S_->Require<CompositeVector,CompositeVectorSpace>(conserved_key_, tag_current_, name_);

  Key molar_dens_key = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
  S_->Require<CompositeVector,CompositeVectorSpace>(molar_dens_key, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(molar_dens_key, tag_next_);

  Key mass_dens_key = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
  S_->Require<CompositeVector,CompositeVectorSpace>(mass_dens_key, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(mass_dens_key, tag_next_);

  Key gas_dens_key = Keys::readKey(*plist_, domain_, "molar density gas", "molar_density_gas");
  S_->Require<CompositeVector,CompositeVectorSpace>(gas_dens_key, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(gas_dens_key, tag_next_);

  // -- thermal conductivity
  // move evaluator from PK plist to State
  if (!S_->HasEvaluator(conductivity_key_, tag_next_)) {
    // only get thermal conductivity if it wasn't set in three-phase_energy
    if (plist_->isSublist("thermal conductivity evaluator")) {
      auto& tcm_plist = S_->GetEvaluatorList(conductivity_key_);
      tcm_plist.setParameters(plist_->sublist("thermal conductivity evaluator"));
      tcm_plist.set("evaluator type", "two-phase thermal conductivity");
    }
    S_->Require<CompositeVector,CompositeVectorSpace>(conductivity_key_, tag_next_).SetMesh(mesh_)
      ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(conductivity_key_, tag_next_);
  }

}

} // namespace
} // namespace
