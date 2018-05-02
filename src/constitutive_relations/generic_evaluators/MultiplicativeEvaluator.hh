/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_MULTIPLICATIVE_EVALUATOR_
#define AMANZI_RELATIONS_MULTIPLICATIVE_EVALUATOR_

#include "factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

class MultiplicativeEvaluator : public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:
  // constructor format for all derived classes
  explicit
  MultiplicativeEvaluator(Teuchos::ParameterList& plist);
  MultiplicativeEvaluator(const MultiplicativeEvaluator& other);

  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondary
  void Evaluate_(const State& S,
                      CompositeVector& result);
  void EvaluatePartialDerivative_(const State& S,
                                       const Key& wrt_key, const Key& wrt_tag, CompositeVector& result);

 protected:
  double coef_;
  
 private:
  static Utils::RegisteredFactory<Evaluator,MultiplicativeEvaluator> factory_;
};

} // namespace
} // namespace

#endif

