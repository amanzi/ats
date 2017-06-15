/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_TWO_PHASE_HH_
#define PKS_ENERGY_TWO_PHASE_HH_

#include "pk_factory.hh"
#include "energy_base.hh"

namespace Amanzi {

// forward declarations
class MPCDiagonalFlowEnergy;
class MPCCoupledFlowEnergy;
namespace Relations { class EOS; }
namespace Energy { namespace EnergyRelations { class IEM; } }

namespace Energy {

class TwoPhase : public EnergyBase {

public:
  TwoPhase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
           Teuchos::ParameterList& FElist,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~TwoPhase() {}

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

protected:
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  //  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S);


 protected:
  // models for evaluating enthalpy
  Teuchos::RCP<Relations::EOS> eos_liquid_;
  Teuchos::RCP<EnergyRelations::IEM> iem_liquid_;

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
