/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for surface temperature
----------------------------------------------------------------------------- */


#ifndef AMANZI_LAKE_SURFACE_TEMPERATURE_EVALUATOR_HH_
#define AMANZI_LAKE_SURFACE_TEMPERATURE_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeSurfaceTemperatureEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LakeSurfaceTemperatureEvaluator(Teuchos::ParameterList& plist);
  LakeSurfaceTemperatureEvaluator(const LakeSurfaceTemperatureEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  Key temperature_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LakeSurfaceTemperatureEvaluator> factory_;

};

} // namespace
} // namespace

#endif
