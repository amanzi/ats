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
                   const Teuchos::RCP<TreeVector>& solution)
  : PK(FElist, plist, S, solution), EnergyBase(FElist, plist, S, solution)
{}


// -------------------------------------------------------------
// Create the physical evaluators for energy and thermal conductivity
// -------------------------------------------------------------
void
TwoPhase::SetupPhysicalEvaluators_()
{
  // -- thermal conductivity
  // This deals with deprecated location for the TC list (in the PK).  Move it
  // to State
  if (plist_->isSublist("thermal conductivity evaluator")) {
    auto& tcm_plist = S_->GetEvaluatorList(conductivity_key_);
    tcm_plist.setParameters(plist_->sublist("thermal conductivity evaluator"));
    tcm_plist.set("evaluator type", "two-phase thermal conductivity");
  }

  EnergyBase::SetupPhysicalEvaluators_();
}

} // namespace Energy
} // namespace Amanzi
