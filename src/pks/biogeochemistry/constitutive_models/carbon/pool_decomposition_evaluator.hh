/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Evaluates carbon pool turnover.

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
  explicit PoolDecompositionEvaluator(Teuchos::ParameterList& plist);
  PoolDecompositionEvaluator(const PoolDecompositionEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void
  EvaluateField_(const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key,
                                               const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key carbon_key_;
  Key decay_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, PoolDecompositionEvaluator> fac_;
};

} // namespace BGCRelations
} // namespace BGC
} // namespace Amanzi

#endif
