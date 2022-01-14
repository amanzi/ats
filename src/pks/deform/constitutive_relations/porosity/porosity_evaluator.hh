/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow according to a Manning approach.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_DEFORMRELATIONS_POROSITY_EVALUATOR_
#define AMANZI_DEFORMRELATIONS_POROSITY_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"



namespace Amanzi {
namespace Deform {
namespace DeformRelations {

class PorosityEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  PorosityEvaluator(Teuchos::ParameterList& cond_plist);
  PorosityEvaluator(const PorosityEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

private:
  static Utils::RegisteredFactory<Evaluator,PorosityEvaluator> factory_;


};

} //namespace
} //namespace
} //namespace

#endif

