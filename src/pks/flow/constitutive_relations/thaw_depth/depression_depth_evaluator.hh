/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The dynamic subgrid model evaluator gets the subgrid parameters and evolve polygons.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_DEPRESSION_DEPTH_EVALUATOR_
#define AMANZI_FLOWRELATIONS_DEPRESSION_DEPTH_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class DepressionDepthEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  DepressionDepthEvaluator(Teuchos::ParameterList& plist);
  DepressionDepthEvaluator(const DepressionDepthEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;
  
protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
  Key delta_init_key_,delta_evolve_key_,sg_entity_key_;

private:

  static Utils::RegisteredFactory<Evaluator,DepressionDepthEvaluator> reg_;  

};
  
} //namespace
} //namespace 

#endif
