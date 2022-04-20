/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "surface_relperm_evaluator.hh"
#include "surface_relperm_model_factory.hh"

namespace Amanzi {
namespace Flow {

SurfaceRelPermEvaluator::SurfaceRelPermEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // create the model
  SurfaceRelPermModelFactory fac;
  model_ = fac.createModel(plist_.sublist("surface rel perm model"));
 
  Key domain;
  if (my_key_.empty()){
    domain = "surface";
    my_key_ = "surface-relative_permeability";
  }
  else
    domain = Keys::getDomain(my_key_);

  // set up the height dependency
  h_key_ = Keys::readKey(plist_, domain, "pressure key", "pressure");

  dependencies_.insert(h_key_);

  // set up the temperature dependency
  is_temp_ = model_->TemperatureDependent();
  if (is_temp_) {
    uf_key_ = Keys::readKey(plist_, domain, "unfrozen fraction", "unfrozen_fraction");
    dependencies_.insert(uf_key_);
  }


}


SurfaceRelPermEvaluator::SurfaceRelPermEvaluator(const SurfaceRelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    is_temp_(other.is_temp_),
    uf_key_(other.uf_key_),
    h_key_(other.h_key_),
    model_(other.model_) {}


Teuchos::RCP<FieldEvaluator>
SurfaceRelPermEvaluator::Clone() const {
  return Teuchos::rcp(new SurfaceRelPermEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void SurfaceRelPermEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  if (is_temp_) {
    Teuchos::RCP<const CompositeVector> uf = S->GetFieldData(uf_key_);
    Teuchos::RCP<const CompositeVector> h = S->GetFieldData(h_key_);

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& uf_v = *uf->ViewComponent(*comp,false);
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->SurfaceRelPerm(uf_v[0][i], h_v[0][i]);
      }
    }

  } else {
    Teuchos::RCP<const CompositeVector> h = S->GetFieldData(h_key_);

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->SurfaceRelPerm(0., h_v[0][i]);
      }
    }
  }
}


void SurfaceRelPermEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(0);
}


} //namespace
} //namespace

