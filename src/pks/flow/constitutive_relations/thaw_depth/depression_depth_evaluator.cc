/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The dynamic subgrid model evaluator gets the subgrid parameters and evolve polygons.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "depression_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
  

DepressionDepthEvaluator::DepressionDepthEvaluator(Teuchos::ParameterList& plist) :
   EvaluatorSecondaryMonotypeCV(plist){

  std::string domain_name=Keys::getDomain(my_key_);
  if (my_key_.empty())
    my_key_ = Keys::getKey(domain_name, "depression_depth");
  delta_init_key_ = plist_.get<std::string>("depression depth initial key");
  dependencies_.insert(delta_init_key_);
  delta_evolve_key_ = plist_.get<std::string>("depression depth evolution key");
  dependencies_.insert(delta_evolve_key_);
  sg_entity_key_ = plist_.get<std::string>("subgrid entity key");
  dependencies_.insert(sg_entity_key_);
}

  
DepressionDepthEvaluator::DepressionDepthEvaluator(const DepressionDepthEvaluator& other) :
  EvaluatorSecondaryMonotypeCV(other),
    delta_init_key_(other.delta_init_key_),
  delta_evolve_key_(other.delta_evolve_key_),
  sg_entity_key_(other.sg_entity_key_) {}


Teuchos::RCP<Evaluator>
DepressionDepthEvaluator::Clone() const {
  return Teuchos::rcp(new DepressionDepthEvaluator(*this));
}
void DepressionDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
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
    
void DepressionDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result){}

void
DepressionDepthEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  AMANZI_ASSERT(my_key_ != std::string(""));
   
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT,  my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  if (my_fac->Mesh() != Teuchos::null) {
    Teuchos::RCP<CompositeVectorSpace> dep_fac =
      Teuchos::rcp(new CompositeVectorSpace());
    dep_fac->SetMesh(my_fac->Mesh());
    dep_fac->AddComponent("cell", AmanziMesh::CELL, 1);
    dep_fac->SetGhosted(true);
    for (const auto& key : dependencies_) { 
      Teuchos::RCP<CompositeVectorSpace> fac = S->Require<CompositeVector,CompositeVectorSpace>(key, Tags::NEXT);
      fac->Update(*dep_fac);
    }
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}
  
} //namespace
} //namespace



  
