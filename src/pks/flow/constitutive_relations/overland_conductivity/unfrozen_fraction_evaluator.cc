/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_fraction_model.hh"
#include "unfrozen_fraction_evaluator.hh"

namespace Amanzi {
namespace Flow {

UnfrozenFractionEvaluator::UnfrozenFractionEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  Key domain = Keys::getDomain(my_key_);

  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain,"temperature"));
  dependencies_.insert(temp_key_);

  if (my_key_ == std::string("")) {
    my_key_ = "surface-unfrozen_fraction";

  }

  // create the model, hard-coded until we have a 2nd model
  AMANZI_ASSERT(plist_.isSublist("unfrozen fraction model"));
  Teuchos::ParameterList sublist = plist_.sublist("unfrozen fraction model");
  model_ = Teuchos::rcp(new UnfrozenFractionModel(sublist));
}


UnfrozenFractionEvaluator::UnfrozenFractionEvaluator(const UnfrozenFractionEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    temp_key_(other.temp_key_),
    model_(other.model_) {}


Teuchos::RCP<Evaluator>
UnfrozenFractionEvaluator::Clone() const {
  return Teuchos::rcp(new UnfrozenFractionEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void UnfrozenFractionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->UnfrozenFraction(temp_v[0][i]);
    }
  }
}


void UnfrozenFractionEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  AMANZI_ASSERT(wrt_key == temp_key_);
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->DUnfrozenFractionDT(temp_v[0][i]);
    }
  }
}


} //namespace
} //namespace

