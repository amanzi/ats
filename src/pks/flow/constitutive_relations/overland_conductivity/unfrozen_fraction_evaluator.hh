/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen fraction model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class UnfrozenFractionModel;

class UnfrozenFractionEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  UnfrozenFractionEvaluator(Teuchos::ParameterList& plist);
  UnfrozenFractionEvaluator(const UnfrozenFractionEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<const UnfrozenFractionModel> get_Model() const { return model_; }
  Teuchos::RCP<UnfrozenFractionModel> get_Model() { return model_; }

protected:
  Teuchos::RCP<UnfrozenFractionModel> model_;
  Key temp_key_;

private:
  static Utils::RegisteredFactory<Evaluator,UnfrozenFractionEvaluator> fac_;


};

} //namespace
} //namespace

#endif

