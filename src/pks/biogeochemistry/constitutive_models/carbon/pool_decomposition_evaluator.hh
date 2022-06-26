/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates carbon pool turnover.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_BGCRELATIONS_POOL_DECOMP_HH_
#define AMANZI_BGCRELATIONS_POOL_DECOMP_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

class PoolDecompositionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit
  PoolDecompositionEvaluator(Teuchos::ParameterList& plist);
  PoolDecompositionEvaluator(const PoolDecompositionEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

protected:
  Key carbon_key_;
  Key decay_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,PoolDecompositionEvaluator> fac_;



};

} // namespace
} // namespace
} // namespace

#endif
