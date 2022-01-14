/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen effective depth, h * eta

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_EFFECTIVE_DEPTH_EVALUATOR_
#define AMANZI_FLOWRELATIONS_UNFROZEN_EFFECTIVE_DEPTH_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class UnfrozenEffectiveDepthModel;

class UnfrozenEffectiveDepthEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  UnfrozenEffectiveDepthEvaluator(Teuchos::ParameterList& plist);
  UnfrozenEffectiveDepthEvaluator(const UnfrozenEffectiveDepthEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

protected:
  Key uf_key_;
  Key depth_key_;
  double alpha_;

 private:
  static Utils::RegisteredFactory<Evaluator,UnfrozenEffectiveDepthEvaluator> fac_;


};

} //namespace
} //namespace

#endif
