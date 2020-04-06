/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! An advection-diffusion equation for energy in two phases.

/*!

This is simply a subsurface energy equation that places a few more requirements
on the base class.  It could probably go away if we refactor to remove
hard-coded evaluators.

.. _energy-two-phase-pk-spec:
.. admonition:: energy-two-phase-pk-spec

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
namespace Relations { class EOS; }

namespace Energy {

class IEM;


class TwoPhase : public EnergyBase {

public:

  TwoPhase(Teuchos::ParameterList& FElist,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~TwoPhase() {}

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

protected:
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  //  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S);


 protected:
  // models for evaluating enthalpy
  Teuchos::RCP<Amanzi::Relations::EOS> eos_liquid_;
  Teuchos::RCP<Energy::IEM> iem_liquid_;

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
