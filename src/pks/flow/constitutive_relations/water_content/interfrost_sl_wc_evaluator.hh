/*
  The interfrost sl water content evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_INTERFROST_SL_WC_EVALUATOR_HH_
#define AMANZI_FLOW_INTERFROST_SL_WC_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostSlWcModel;

class InterfrostSlWcEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  InterfrostSlWcEvaluator(Teuchos::ParameterList& plist);
  InterfrostSlWcEvaluator(const InterfrostSlWcEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<InterfrostSlWcModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

 protected:
  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key ni_key_;
  Key cv_key_;

  Teuchos::RCP<InterfrostSlWcModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,InterfrostSlWcEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
