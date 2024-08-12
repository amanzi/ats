/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

Evaluator for enthalpy.
----------------------------------------------------------------------------- */


#ifndef AMANZI_LAKE_ENTHALPY_EVALUATOR_HH_
#define AMANZI_LAKE_ENTHALPY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace LakeThermo {

class LakeEnthalpyEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  LakeEnthalpyEvaluator(Teuchos::ParameterList& plist);
  LakeEnthalpyEvaluator(const LakeEnthalpyEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:

  Key pres_key_;
  Key dens_key_;
  Key ie_key_;
  bool include_work_;

 private:
  static Utils::RegisteredFactory<Evaluator,LakeEnthalpyEvaluator> factory_;

};

} // namespace
} // namespace

#endif
