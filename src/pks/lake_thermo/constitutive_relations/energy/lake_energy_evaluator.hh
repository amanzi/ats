/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

Evaluator for energy, e = cp*rho*T.
----------------------------------------------------------------------------- */


#ifndef AMANZI_LAKE_ENERGY_EVALUATOR_HH_
#define AMANZI_LAKE_ENERGY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeEnergyEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  LakeEnergyEvaluator(Teuchos::ParameterList& plist);
  LakeEnergyEvaluator(const LakeEnergyEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:

  Key temperature_key_;
  Key density_key_;
  Key heat_capacity_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,LakeEnergyEvaluator> factory_;

};

} // namespace
} // namespace

#endif
