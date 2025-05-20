/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Modifies the energy equation to provide two-phase energy (water + vapor/air).

`"PK type`" = `"two-phase energy`"

.. _pk-two-phase-energy-spec:
.. admonition:: pk-two-phase-energy-spec

   INCLUDES:

   - ``[energy-pk-spec]``  See `Energy Base PK`_

*/

#ifndef PKS_ENERGY_TWO_PHASE_HH_
#define PKS_ENERGY_TWO_PHASE_HH_

#include "PK_Factory.hh"
#include "energy_base.hh"

namespace Amanzi {

// forward declarations
class MPCDiagonalFlowEnergy;
class MPCCoupledFlowEnergy;

namespace Energy {

class TwoPhase : public EnergyBase {
 public:
  TwoPhase(Teuchos::ParameterList& FElist,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~TwoPhase() override {}

 protected:
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_() override;

 private:
  // factory registration
  static RegisteredPKFactory<TwoPhase> reg_;

  // Energy has a friend in couplers...
  friend class Amanzi::MPCCoupledFlowEnergy;
  friend class Amanzi::MPCDiagonalFlowEnergy;
};

} // namespace Energy
} // namespace Amanzi

#endif
