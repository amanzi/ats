/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow subgrid model.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_SUBGRID_EVALUATOR_
#define AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_SUBGRID_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivitySubgridEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  OverlandConductivitySubgridEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivitySubgridEvaluator(const OverlandConductivitySubgridEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<ManningConductivityModel> get_Model() { return model_; }

private:
  Teuchos::RCP<ManningConductivityModel> model_;

  Key slope_key_;
  Key coef_key_;
  Key dens_key_;
  Key depth_key_;
  Key drag_exp_key_;
  Key frac_cond_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandConductivitySubgridEvaluator> factory_;
};

} //namespace
} //namespace

#endif

