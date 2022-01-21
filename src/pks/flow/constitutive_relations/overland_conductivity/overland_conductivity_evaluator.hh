/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_EVALUATOR_
#define AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivityEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  OverlandConductivityEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;

  Teuchos::RCP<ManningConductivityModel> get_Model() { return model_; }

private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key depth_key_;
  Key slope_key_;
  Key coef_key_;
  Key dens_key_;
  double dt_swe_factor_;
  bool dens_;

 private:
  static Utils::RegisteredFactory<Evaluator,OverlandConductivityEvaluator> factory_;
};

} //namespace
} //namespace

#endif

