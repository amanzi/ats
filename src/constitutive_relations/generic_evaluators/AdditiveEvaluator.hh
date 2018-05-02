/*
  AdditiveEvaluator is the generic evaluator for adding N other fields.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_ADDITIVE_EVALUATOR_
#define AMANZI_RELATIONS_ADDITIVE_EVALUATOR_

#include "factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

  class AdditiveEvaluator : public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:
  // constructor format for all derived classes
  explicit
  AdditiveEvaluator(Teuchos::ParameterList& plist);

  AdditiveEvaluator(const AdditiveEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondary
  void Evaluate_(const State& S,
                      CompositeVector& result);
  void EvaluatePartialDerivative_(const State& S,
                                       const Key& wrt_key, const Key& wrt_tag, CompositeVector& result);

 protected:
  std::map<Key, double> coefs_;

 private:
  static Utils::RegisteredFactory<Evaluator,AdditiveEvaluator> factory_;
};

} // namespace
} // namespace

#endif

