/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a water/ice heat capacity model

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#ifndef AMANZI_LAKE_HC_EVALUATOR_HH_
#define AMANZI_LAKE_HC_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeHeatCapacityEvaluator :
    public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  LakeHeatCapacityEvaluator(Teuchos::ParameterList& plist);
  LakeHeatCapacityEvaluator(const LakeHeatCapacityEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldModel
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  // dependencies

  double cw, ci;

  Key temperature_key_;
  Key water_content_key_;
  Key ice_content_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,LakeHeatCapacityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
