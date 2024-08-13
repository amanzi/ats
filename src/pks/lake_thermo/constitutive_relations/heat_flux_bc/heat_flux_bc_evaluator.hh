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
#include "EvaluatorSecondaryMonotype.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace LakeThermo {

class HeatFluxBCEvaluator :
    public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  HeatFluxBCEvaluator(Teuchos::ParameterList& plist);
  HeatFluxBCEvaluator(const HeatFluxBCEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldModel
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

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
  static Utils::RegisteredFactory<Evaluator,HeatFluxBCEvaluator> factory_;

};

} // namespace
} // namespace

#endif
