/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow subgrid model.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_SUBGRID_EVALUATOR_
#define AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_SUBGRID_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivitySubgridEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  OverlandConductivitySubgridEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivitySubgridEvaluator(const OverlandConductivitySubgridEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<ManningConductivityModel> get_Model() { return model_; }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key slope_key_;
  Key coef_key_;
  Key dens_key_;
  Key mobile_depth_key_;
  Key drag_exp_key_;
  Key frac_cond_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,OverlandConductivitySubgridEvaluator> factory_;
};

} //namespace
} //namespace

#endif

