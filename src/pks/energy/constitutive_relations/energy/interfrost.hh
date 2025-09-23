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

#ifndef PKS_ENERGY_THREE_PHASE_HH_
#define PKS_ENERGY_THREE_PHASE_HH_

#include "pk_factory_ats.hh"
#include "three_phase.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class Interfrost : public ThreePhase {
 public:
  Interfrost(Teuchos::ParameterList& FElist,
             const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& solution)
    : //PKDefaultBase(plist, FElist, solution),
      ThreePhase(FElist, plist, S, solution)
  {
    plist_ = plist;
    solution_ = solution;
  }


  // Virtual destructor
  virtual ~Interfrost() {}

 protected:
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& f);
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

 private:
  static RegisteredPKFactory_ATS<Interfrost> reg_;

  friend class MPCCoupledFlowEnergy;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
