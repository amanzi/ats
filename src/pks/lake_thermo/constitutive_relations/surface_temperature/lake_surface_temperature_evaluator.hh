/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for surface temperature
----------------------------------------------------------------------------- */


#ifndef AMANZI_LAKE_SURFACE_TEMPERATURE_EVALUATOR_HH_
#define AMANZI_LAKE_SURFACE_TEMPERATURE_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeSurfaceTemperatureEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  LakeSurfaceTemperatureEvaluator(Teuchos::ParameterList& plist);
  LakeSurfaceTemperatureEvaluator(const LakeSurfaceTemperatureEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:

  Key temperature_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,LakeSurfaceTemperatureEvaluator> factory_;

};

} // namespace
} // namespace

#endif
