/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_EFFECTIVE_PRESSURE_EVALUATOR_HH_
#define AMANZI_EFFECTIVE_PRESSURE_EVALUATOR_HH_

#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

class EffectivePressureEvaluator : public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:

  // constructor format for all derived classes
  explicit
  EffectivePressureEvaluator(Teuchos::ParameterList& ep_plist);

  EffectivePressureEvaluator(const EffectivePressureEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from SecondaryVariableEvaluator
  virtual void Evaluate_(const State& S,
          CompositeVector& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Key& wrt_tag, CompositeVector& result) override;

 protected:

  // PList
  Teuchos::ParameterList ep_plist_;

  // Keys for fields
  // dependencies
  Key pres_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,EffectivePressureEvaluator> factory_;
};

} // namespace
} // namespace

#endif
