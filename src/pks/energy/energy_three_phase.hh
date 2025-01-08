/*
  Copyright 2010-202x held jointly by participating institutions.
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

.. _energy_three_phase_pk-spec:
.. admonition:: energy_three_phase_pk-spec

    INCLUDES:

    - ``[energy_two_phase_pk-spec]`` See  `Two-Phase subsurface Energy PK`_

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
             const Teuchos::RCP<TreeVector>& solution);

  virtual void Initialize() override;

 protected:
  virtual void SetupPhysicalEvaluators_() override;

 private:
  static RegisteredPKFactory<ThreePhase> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
