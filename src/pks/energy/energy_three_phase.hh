/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! An advection-diffusion equation for energy in three phases.

/*!

This is simply a subsurface energy equation that places a few more requirements
on the base class.  It could probably go away if we refactor to remove
hard-coded evaluators.

.. _energy-three-phase-pk-spec:
.. admonition:: energy-three-phase-pk-spec

    INCLUDES:

    - ``[energy-two-phase-pk-spec]`` See  `Two-Phase subsurface Energy PK`_

*/

#ifndef PKS_ENERGY_THREE_PHASE_HH_
#define PKS_ENERGY_THREE_PHASE_HH_

#include "energy_two_phase.hh"

namespace Amanzi {
namespace Energy {

class ThreePhase : public TwoPhase {

public:

  ThreePhase(Teuchos::ParameterList& FElist,
             const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    TwoPhase(FElist, plist, S, solution) {}

  // Virtual destructor
  virtual ~ThreePhase() {}

  virtual void Initialize(const Teuchos::Ptr<State>& S);
  
protected:
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

private:
  static RegisteredPKFactory<ThreePhase> reg_;

  friend  class MPCCoupledFlowEnergy;
};

} // namespace Energy
} // namespace Amanzi

#endif
