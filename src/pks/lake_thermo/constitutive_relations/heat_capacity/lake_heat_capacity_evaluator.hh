/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a water/ice heat capacity model

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#ifndef AMANZI_LAKE_HC_EVALUATOR_HH_
#define AMANZI_LAKE_HC_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeHeatCapacityEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  LakeHeatCapacityEvaluator(Teuchos::ParameterList& plist);
  LakeHeatCapacityEvaluator(const LakeHeatCapacityEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // dependencies

  double cw, ci;

  Key temperature_key_;
  Key water_content_key_;
  Key ice_content_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LakeHeatCapacityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
