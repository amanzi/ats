/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes the maximum thaw depth 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "max_thaw_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
  

MaxThawDepthEvaluator::MaxThawDepthEvaluator(Teuchos::ParameterList& plist) :
   SecondaryVariableFieldEvaluator(plist){
  threshold_td_ = plist_.get<double>("threshold value", 0.4);

  std::string domain_name=Keys::getDomain(my_key_);
  my_key_ = Keys::getKey(domain_name, "maximum_thaw_depth");
  td_key_ = plist_.get<std::string>("thaw depth key", Keys::getKey(domain_name, "thaw_depth"));
  dependencies_.insert(td_key_);
}

  
MaxThawDepthEvaluator::MaxThawDepthEvaluator(const MaxThawDepthEvaluator& other) :
  SecondaryVariableFieldEvaluator(other),
    td_key_(other.td_key_),
    threshold_td_(other.threshold_td_){}

Teuchos::RCP<FieldEvaluator>
MaxThawDepthEvaluator::Clone() const {
  return Teuchos::rcp(new MaxThawDepthEvaluator(*this));
}
void MaxThawDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& thawdepth_c = *S->GetFieldData(td_key_)->ViewComponent("cell", false);
  Key domain = Keys::getDomain(td_key_);
  assert(!domain.empty());

  int ncells = res_c.MyLength();
  AMANZI_ASSERT(ncells ==1);
  for (int c=0; c!=ncells; c++){
    res_c[0][c] = std::max<double>(thawdepth_c[0][c], res_c[0][c]);
  }

}
    
void MaxThawDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result){}

} //namespace
} //namespace



  
