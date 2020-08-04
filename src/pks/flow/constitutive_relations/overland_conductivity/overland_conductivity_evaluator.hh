/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_EVALUATOR_
#define AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel;

class OverlandConductivityEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  OverlandConductivityEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

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
  static Utils::RegisteredFactory<FieldEvaluator,OverlandConductivityEvaluator> factory_;
};

} //namespace
} //namespace

#endif

