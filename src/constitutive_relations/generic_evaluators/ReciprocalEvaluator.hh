/*
  ReciprocalEvaluator is the generic evaluator for dividing two vectors.

  Authors: Daniil Svyatsky
*/

#ifndef AMANZI_RELATIONS_RECIPROCAL_EVALUATOR_
#define AMANZI_RELATIONS_RECIPROCAL_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class ReciprocalEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  explicit
  ReciprocalEvaluator(Teuchos::ParameterList& plist);
  ReciprocalEvaluator(const ReciprocalEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  void Evaluate_(const State& S,
                 const std::vector<CompositeVector*>& result) override;
  void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

 protected:
  double coef_;
  bool positive_;
  Key reciprocal_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,ReciprocalEvaluator> factory_;
};

} // namespace
} // namespace

#endif

