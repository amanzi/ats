/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for enthalpy.
----------------------------------------------------------------------------- */


#include "enthalpy_evaluator.hh"

namespace Amanzi {
namespace Energy {

EnthalpyEvaluator::EnthalpyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_.empty()) {

    my_key_ = plist_.get<std::string>("enthalpy key", "surface-enthalpy_liquid");
  }

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);
  include_work_ = plist_.get<bool>("include work term", false);

  // -- pressure
  if (include_work_) {
    pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
    dependencies_.insert(pres_key_);

    dens_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
    dependencies_.insert(dens_key_);
  }

  ie_key_ = Keys::readKey(plist_, domain_name, "internal energy liquid", "internal_energy_liquid");
  dependencies_.insert(ie_key_);

};

EnthalpyEvaluator::EnthalpyEvaluator(const EnthalpyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    dens_key_(other.dens_key_),
    ie_key_(other.ie_key_),
    include_work_(other.include_work_) {};

Teuchos::RCP<FieldEvaluator>
EnthalpyEvaluator::Clone() const {
  return Teuchos::rcp(new EnthalpyEvaluator(*this));
};


void EnthalpyEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::OSTab tab = vo_->getOSTab();
  Teuchos::RCP<const CompositeVector> u_l = S->GetPtr<CompositeVector>(ie_key_);
  *result = *u_l;


  if (include_work_) {
    Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(pres_key_);
    Teuchos::RCP<const CompositeVector> n_l = S->GetPtr<CompositeVector>(dens_key_);

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp,false);
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] += 1.e-6*pres_v[0][i]/nl_v[0][i];
      }
    }
  }
};


void EnthalpyEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // not implemented
  if (wrt_key == ie_key_) {
    result->PutScalar(1.);
  } else if (wrt_key ==pres_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> n_l = S->GetPtr<CompositeVector>(dens_key_);

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] = 1.e-6/nl_v[0][i];
      }
    }

  } else if (wrt_key ==dens_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(pres_key_);
    Teuchos::RCP<const CompositeVector> n_l = S->GetPtr<CompositeVector>(dens_key_);

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp,false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] = -1.e-6*pres_v[0][i]/std::pow(nl_v[0][i], 2);
      }
    }
  }
};


} //namespace
} //namespace
