/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_LAKE_TC_EVALUATOR_HH_
#define AMANZI_LAKE_TC_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace LakeThermo {

class ThermalConductivityEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  ThermalConductivityEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityEvaluator(const ThermalConductivityEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // dependencies

  bool ice_cover_ = false;

  double V_wind_;
  double V_wind_0_;
  double K_max_;
  double K_0_;

  Key temperature_key_;
  Key depth_key_;


//  Key uf_key_;
//  Key height_key_;

//  double K_liq_;
//  double K_ice_;
//  double min_K_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,ThermalConductivityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
