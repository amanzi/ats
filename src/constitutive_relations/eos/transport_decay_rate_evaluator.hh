/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

/*
  The elevation evaluator calculates temperature-dependent decay rate.

*/

#ifndef AMANZI_FLOWRELATIONS_TRANSPORTDECAY_EVALUATOR_
#define AMANZI_FLOWRELATIONS_TRANSPORTDECAY_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class TransportDecayRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit TransportDecayRateEvaluator(Teuchos::ParameterList& plist);
  TransportDecayRateEvaluator(const TransportDecayRateEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override
  {}
  
  double Func_Temp(double temp, double temp_ref, double q10) const;

 protected:
  Key temp_key_;
  Key domain_;

  double q10_, temp_ref_, decay_ref_;

 private:
  static Utils::RegisteredFactory<Evaluator, TransportDecayRateEvaluator> reg_;
};

} // namespace Relations
} // namespace Amanzi

#endif
