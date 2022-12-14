/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_WATER_CONTENT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_WATER_CONTENT_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class OverlandPressureWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit OverlandPressureWaterContentEvaluator(Teuchos::ParameterList& plist);
  OverlandPressureWaterContentEvaluator(const OverlandPressureWaterContentEvaluator& other) =
    default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  void InitializeFromPlist_();

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key pres_key_, cv_key_;

  double M_;
  bool bar_; // bar'd variable indicates this is potentially negative for
             // pressures less than atmospheric
  double rollover_;

 private:
  static Utils::RegisteredFactory<Evaluator, OverlandPressureWaterContentEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
