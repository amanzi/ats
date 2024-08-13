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
#include "EvaluatorSecondaryMonotype.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace LakeThermo {

class ThermalConductivityEvaluator :
    public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  ThermalConductivityEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityEvaluator(const ThermalConductivityEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldModel
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

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
  static Utils::RegisteredFactory<Evaluator,ThermalConductivityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
