/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a water/ice heat capacity model

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#ifndef AMANZI_LAKE_ER_EVALUATOR_HH_
#define AMANZI_LAKE_ER_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeEvaporationRateEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  LakeEvaporationRateEvaluator(Teuchos::ParameterList& plist);
  LakeEvaporationRateEvaluator(const LakeEvaporationRateEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // dependencies

  Key temperature_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LakeEvaporationRateEvaluator> factory_;

};

} // namespace
} // namespace

#endif
