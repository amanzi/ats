/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_
#define AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "thermal_conductivity_twophase.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class ThermalConductivityTwoPhaseEvaluator :
    public EvaluatorSecondaryMonotypeCV {

 public:
  using RegionModelPair = std::pair<std::string,Teuchos::RCP<ThermalConductivityTwoPhase> >;

  // constructor format for all derived classes
  ThermalConductivityTwoPhaseEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityTwoPhaseEvaluator(const ThermalConductivityTwoPhaseEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  std::vector<RegionModelPair> tcs_;

  // Keys for fields
  // dependencies
  Key poro_key_;
  Key sat_key_;
};

} // namespace
} // namespace

#endif
