/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon  (coonet@ornl.gov)
*/

/*
  An evaluator for compute the divergence of a flux field.
*/
#pragma once

#include "dbc.hh"
#include "Evaluator_Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class FluxDivergenceEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit FluxDivergenceEvaluator(Teuchos::ParameterList& plist);
  FluxDivergenceEvaluator(const FluxDivergenceEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // turn off all derivatives manually
  virtual bool IsDifferentiableWRT(const State& S,
          const Key& wrt_key,
          const Tag& wrt_tag) const override {
    return false;
  }

protected:
  // custom ensure compatibility as all data is not just on the same components
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key,
          const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override {
    AMANZI_ASSERT(false);
  }

 protected:
  Key flux_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, FluxDivergenceEvaluator> reg_;
};


} // namespace Relations
} // namespace Amanzi


