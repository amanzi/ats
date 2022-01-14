/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "effective_height_model.hh"
#include "effective_height_evaluator.hh"


namespace Amanzi {
namespace Flow {

EffectiveHeightEvaluator::EffectiveHeightEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {
  // my keys are for saturation and rel perm.
  if (my_key_ == "")
    my_key_ = plist_.get<std::string>("effective height key", "effective_height");

  // my dependencies
  height_key_ = plist_.get<std::string>("height key", "ponded_depth");
  dependencies_.insert(height_key_);

  // model
  Teuchos::ParameterList model_plist = plist_.sublist("effective height model parameters");
  model_ = Teuchos::rcp(new EffectiveHeightModel(model_plist));
}


EffectiveHeightEvaluator::EffectiveHeightEvaluator(const EffectiveHeightEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    height_key_(other.height_key_),
    model_(other.model_) {}


Teuchos::RCP<Evaluator>
EffectiveHeightEvaluator::Clone() const {
  return Teuchos::rcp(new EffectiveHeightEvaluator(*this));
}

void EffectiveHeightEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> height = S->GetPtr<CompositeVector>(height_key_);

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& height_v = *(height->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = model_->EffectiveHeight(height_v[0][id]);
    }
  }
}


void EffectiveHeightEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(wrt_key == height_key_);

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> height = S->GetPtr<CompositeVector>(height_key_);

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& height_v = *(height->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = model_->DEffectiveHeightDHeight(height_v[0][id]);
    }
  }

}



} //namespace
} //namespace
