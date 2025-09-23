/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

Modifies the energy equation to provide three-phase energy (water + ice + vapor/air).

`"PK type`" = `"three-phase energy`"

.. _pk-three-phase-energy-spec:
.. admonition:: pk-three-phase-energy-spec

   INCLUDES:

   - ``[pk-two-phase-energy-spec]`` See  :ref:`Two-Phase Energy PK`

*/

#ifndef PKS_ENERGY_THREE_PHASE_HH_
#define PKS_ENERGY_THREE_PHASE_HH_

#include "energy_two_phase.hh"

namespace Amanzi {
namespace ATS_Physics {    
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
} // namespace ATS_Physics        
} // namespace Amanzi

#endif
