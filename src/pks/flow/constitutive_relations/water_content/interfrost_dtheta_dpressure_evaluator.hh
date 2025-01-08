/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The interfrost dtheta_dpressure evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Interfrost water content portion sl.

*/

#ifndef AMANZI_FLOW_INTERFROST_DTHETA_DPRESSURE_EVALUATOR_HH_
#define AMANZI_FLOW_INTERFROST_DTHETA_DPRESSURE_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostDthetaDpressureModel;

class InterfrostDthetaDpressureEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit InterfrostDthetaDpressureEvaluator(Teuchos::ParameterList& plist);
  InterfrostDthetaDpressureEvaluator(const InterfrostDthetaDpressureEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  Teuchos::RCP<InterfrostDthetaDpressureModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key nl_key_;
  Key sl_key_;
  Key phi_key_;

  Teuchos::RCP<InterfrostDthetaDpressureModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, InterfrostDthetaDpressureEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
