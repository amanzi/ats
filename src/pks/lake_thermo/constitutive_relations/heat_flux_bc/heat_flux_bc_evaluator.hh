/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a surface heat flux in lake model

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#ifndef AMANZI_LAKE_HEAT_FLUX_BC_EVALUATOR_HH_
#define AMANZI_LAKE_HEAT_FLUX_BC_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace LakeThermo {

class HeatFluxBCEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  HeatFluxBCEvaluator(Teuchos::ParameterList& plist);
  HeatFluxBCEvaluator(const HeatFluxBCEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // dependencies

  bool ice_cover_ = false;

  double SS;
  double alpha_w;
  double alpha_i;
  double E_a;
  double E_s;
  double H;
  double LE;

  Key temperature_key_;
  Key conductivity_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,HeatFluxBCEvaluator> factory_;

};

} // namespace
} // namespace

#endif
