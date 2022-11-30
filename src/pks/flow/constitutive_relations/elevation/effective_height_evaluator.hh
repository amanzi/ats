/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining effective_height(height), which is a smoothing
  term near 0 height.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Computes ponded depth from surface water pressure using a smoothed term to make
derivative smooth near 0.  This is pretty much never used anymore.

`"evaluator type`" = `"effective height`"

.. _effective-height-evaluator-spec:
.. admonition:: effective-height-evaluator-spec

   * `"smoothing width [m]`" ``[double]`` **0.01** the width over which smoothing
     is applied.

   KEYS:

   - `"height`" The unsmoothed ponded depth
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
  explicit EffectiveHeightEvaluator(Teuchos::ParameterList& plist);
  EffectiveHeightEvaluator(const EffectiveHeightEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<EffectiveHeightModel> get_Model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key height_key_;

  Teuchos::RCP<EffectiveHeightModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, EffectiveHeightEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
