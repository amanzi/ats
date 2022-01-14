/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The dynamic subgrid model evaluator gets the subgrid parameters and evolve polygons.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "excluded_volume_evaluator.hh"

namespace Amanzi {
namespace Flow {
  

ExcludedVolumeEvaluator::ExcludedVolumeEvaluator(Teuchos::ParameterList& plist) :
   EvaluatorSecondaryMonotypeCV(plist){

  std::string domain_name=Keys::getDomain(my_key_);
  if (my_key_.empty())
    my_key_ = Keys::getKey(domain_name, "excluded_volume");
  delta_init_key_ = plist_.get<std::string>("excluded volume initial key");
  dependencies_.insert(delta_init_key_);
  delta_evolve_key_ = plist_.get<std::string>("excluded volume evolution key");
  dependencies_.insert(delta_evolve_key_);
  sg_entity_key_ = plist_.get<std::string>("subgrid entity key");
  dependencies_.insert(sg_entity_key_);
}

  
ExcludedVolumeEvaluator::ExcludedVolumeEvaluator(const ExcludedVolumeEvaluator& other) :
  EvaluatorSecondaryMonotypeCV(other),
    delta_init_key_(other.delta_init_key_),
  delta_evolve_key_(other.delta_evolve_key_),
  sg_entity_key_(other.sg_entity_key_) {}


Teuchos::RCP<Evaluator>
ExcludedVolumeEvaluator::Clone() const {
  return Teuchos::rcp(new ExcludedVolumeEvaluator(*this));
}
void ExcludedVolumeEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& delta_init_c = *S->Get<CompositeVector>(delta_init_key_).ViewComponent("cell", false);
  const Epetra_MultiVector& delta_evolve_c = *S->Get<CompositeVector>(delta_evolve_key_).ViewComponent("cell",false);
  const Epetra_MultiVector& sg_entity_c = *S->Get<CompositeVector>(sg_entity_key_).ViewComponent("cell",false);

  Key domain = Keys::getDomain(delta_init_key_);
  assert(!domain.empty());

  int ncells = res_c.MyLength();
  for (int c=0; c!=ncells; c++){
    if (sg_entity_c[0][c] == 1) // 1 is for HCP, 0 is for LCP
      res_c[0][c] = delta_init_c[0][c];
    else {
      res_c[0][c] = delta_evolve_c[0][c];
    }
  }

}
    
void ExcludedVolumeEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result){}

} //namespace
} //namespace



  
