/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

Evaluator for water density.
----------------------------------------------------------------------------- */


#ifndef AMANZI_LAKE_DENSITY_EVALUATOR_HH_
#define AMANZI_LAKE_DENSITY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace LakeThermo {

class DensityEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit DensityEvaluator(Teuchos::ParameterList& plist);
  DensityEvaluator(const DensityEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  
  
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:

  Key temperature_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,DensityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
