/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates bioturbation of carbon -- simple diffusion model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_BGCRELATIONS_BIOTURBATION_HH_
#define AMANZI_BGCRELATIONS_BIOTURBATION_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

class BioturbationEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit BioturbationEvaluator(Teuchos::ParameterList& plist);
  BioturbationEvaluator(const BioturbationEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

protected:
  Key carbon_key_;
  Key diffusivity_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,BioturbationEvaluator> fac_;



};

} // namespace
} // namespace
} // namespace

#endif
