/*
  ReciprocalEvaluator is the generic evaluator for dividing two vectors.

  Authors: Daniil Svyatsky
*/

#ifndef AMANZI_RELATIONS_RECIPROCAL_EVALUATOR_
#define AMANZI_RELATIONS_RECIPROCAL_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class ReciprocalEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  ReciprocalEvaluator(Teuchos::ParameterList& plist);
  ReciprocalEvaluator(const ReciprocalEvaluator& other) = default;

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  void EvaluateField_(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<CompositeVector>& result);
  void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  double coef_;
  bool positive_;
  Key reciprocal_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,ReciprocalEvaluator> factory_;
};

} // namespace
} // namespace

#endif

