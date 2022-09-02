/*
  The interfrost denergy_dtemperature evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_INTERFROST_DENERGY_DTEMPERATURE_EVALUATOR_HH_
#define AMANZI_FLOW_INTERFROST_DENERGY_DTEMPERATURE_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostDenergyDtemperatureModel;

class InterfrostDenergyDtemperatureEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  InterfrostDenergyDtemperatureEvaluator(Teuchos::ParameterList& plist);
  InterfrostDenergyDtemperatureEvaluator(const InterfrostDenergyDtemperatureEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<InterfrostDenergyDtemperatureModel> get_model() { return model_; }

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
  Key si_key_;
  Key ni_key_;
  Key rhos_key_;
  Key T_key_;

  Teuchos::RCP<InterfrostDenergyDtemperatureModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,InterfrostDenergyDtemperatureEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
