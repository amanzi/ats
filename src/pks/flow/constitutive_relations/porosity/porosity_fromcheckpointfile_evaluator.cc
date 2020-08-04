/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Reads porosity from checkpoint file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "porosity_fromcheckpointfile_evaluator.hh"


namespace Amanzi {
namespace Flow {

PorosityFromCheckpointFileEvaluator::PorosityFromCheckpointFileEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  std::string domain_name=Keys::getDomain(my_key_);
  
  //poro_key_ = plist_.get<std::string>("base porosity key", Keys::getKey(domain_name,"base_porosity"));
  //dependencies_.insert(poro_key_);

  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("porosity key",
                                      Keys::getKey(domain_name, "porosity"));
  }

  
}


PorosityFromCheckpointFileEvaluator::PorosityFromCheckpointFileEvaluator(const PorosityFromCheckpointFileEvaluator& other) :
  SecondaryVariableFieldEvaluator(other){}
    //poro_key_(other.poro_key_){}

Teuchos::RCP<FieldEvaluator>
PorosityFromCheckpointFileEvaluator::Clone() const {
  return Teuchos::rcp(new PorosityFromCheckpointFileEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void PorosityFromCheckpointFileEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  //  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);

  // evaluate the model
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      std::cout<<"PorosityFrom: "<<"\n";
      result_v[0][id] = result_v[0][id];
    }
  }
}


void PorosityFromCheckpointFileEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  for (CompositeVector::name_iterator comp=result->begin();                                                                                       comp!=result->end(); ++comp) {
    int count = result->size(*comp);
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = 0;
    }
  }
}


} //namespace
} //namespace
