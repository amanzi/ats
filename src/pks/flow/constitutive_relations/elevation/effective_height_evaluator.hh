/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining effective_height(height), which is a smoothing
  term near 0 height.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_EFFECTIVE_HEIGHT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_EFFECTIVE_HEIGHT_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class EffectiveHeightModel;

class EffectiveHeightEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  explicit
  EffectiveHeightEvaluator(Teuchos::ParameterList& plist);
  EffectiveHeightEvaluator(const EffectiveHeightEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;

  Teuchos::RCP<EffectiveHeightModel> get_Model() { return model_; }

 protected:

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key height_key_;

  Teuchos::RCP<EffectiveHeightModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,EffectiveHeightEvaluator> factory_;

};

} //namespace
} //namespace

#endif
