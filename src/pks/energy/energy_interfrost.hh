/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Fully three-phase (air, water, ice) permafrost energy equation, with only
mobile water.

Inherits TwoPhase instead of EnergyBase to pick up the enthalpy from TwoPhase.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_INTERFROST_ENERGY_HH_
#define PKS_ENERGY_INTERFROST_ENERGY_HH_

#include "energy_three_phase.hh"

namespace Amanzi {
namespace ATS_Physics {    
namespace Energy {

class InterfrostEnergy : public ThreePhase {
 public:
  InterfrostEnergy(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution), ThreePhase(FElist, plist, S, solution)
  {}

  // Virtual destructor
  virtual ~InterfrostEnergy() {}

  // -- accumulation term
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

 protected:
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& f) override;
  virtual void SetupPhysicalEvaluators_() override;

 private:
  static RegisteredPKFactory<InterfrostEnergy> reg_;

  friend class MPCCoupledFlowEnergy;
};

} // namespace Energy
} // namespace ATS_Physics          
} // namespace Amanzi

#endif
