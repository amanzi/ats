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

  // -- thermal conductivity
  // move evaluator from PK plist to State
  if (plist_->isSublist("thermal conductivity evaluator")) {
    auto& tcm_plist = S_->GetEvaluatorList(conductivity_key_);
    tcm_plist.setParameters(plist_->sublist("thermal conductivity evaluator"));
    tcm_plist.set("evaluator type", "two-phase thermal conductivity");
  }
  S_->Require<CompositeVector,CompositeVectorSpace>(conductivity_key_, tag_next_).SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(conductivity_key_, tag_next_);

}

} // namespace
} // namespace
