/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "pres_elev_evaluator.hh"

namespace Amanzi {
namespace Flow {

PresElevEvaluator::PresElevEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {
  Key domain = Keys::getDomain(my_key_);

  pres_key_ = plist_.get<std::string>("ponded depth key", Keys::getKey(domain,"ponded_depth"));
  dependencies_.insert(pres_key_);
  elev_key_ = plist_.get<std::string>("elevation key", Keys::getKey(domain,"elevation"));
  dependencies_.insert(elev_key_);
}


PresElevEvaluator::PresElevEvaluator(const PresElevEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    elev_key_(other.elev_key_),
    pres_key_(other.pres_key_) {};

Teuchos::RCP<Evaluator>
PresElevEvaluator::Clone() const {
  return Teuchos::rcp(new PresElevEvaluator(*this));
}


void PresElevEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  // update pressure + elevation
  Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(pres_key_);
  Teuchos::RCP<const CompositeVector> elev = S->GetPtr<CompositeVector>(elev_key_);

  result->Update(1.0, *elev, 1.0, *pres, 0.0);
}


// This is hopefully never called?
void PresElevEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  result->PutScalar(1.0);
}

} //namespace
} //namespace
