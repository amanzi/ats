/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.


  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_SETTLEMENTRATE_EVALUATOR_
#define AMANZI_SETTLEMENTRATE_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class SettlementRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SettlementRateEvaluator(Teuchos::ParameterList& plist);
  SettlementRateEvaluator(const SettlementRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  double tau_d_;
  double ws_;
  double gamma_;
  double lambda_, umax_, xi_;
  double sediment_density_;
  double Cf_;
  Key velocity_key_, sediment_key_;

  static Utils::RegisteredFactory<Evaluator, SettlementRateEvaluator> factory_;
};

} // namespace Amanzi

#endif
