/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Thermal conductivity based on a three-phase (air,liquid,ice) composition of the porous media.

#pragma once

#include "EvaluatorSecondaryMonotype.hh"
#include "thermal_conductivity_threephase.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class ThermalConductivityThreePhaseEvaluator :
    public EvaluatorSecondaryMonotypeCV {

 public:

  typedef std::pair<std::string,Teuchos::RCP<ThermalConductivityThreePhase> > RegionModelPair;

  // constructor format for all derived classes
  ThermalConductivityThreePhaseEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityThreePhaseEvaluator(const ThermalConductivityThreePhaseEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from SecondaryVariableFieldModel
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

 protected:
  std::vector<RegionModelPair> tcs_;

  // Keys for fields
  // dependencies
  Key poro_key_;
  Key sat_key_;
  Key sat2_key_;
  Key temp_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,ThermalConductivityThreePhaseEvaluator> factory_;
};

} // namespace
} // namespace

