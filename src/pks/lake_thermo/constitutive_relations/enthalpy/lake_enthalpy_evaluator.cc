/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for enthalpy.
----------------------------------------------------------------------------- */


#include "lake_enthalpy_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

LakeEnthalpyEvaluator::LakeEnthalpyEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  std::cout << "check in LakeEnthalpyEvaluator" << std::endl;

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  include_work_ = plist_.get<bool>("include work term", true);

  // -- pressure
//  if (include_work_) {
//    pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
//    dependencies_.insert(pres_key_);
//
//    dens_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
//    dependencies_.insert(dens_key_);
//  }

//  ie_key_ = Keys::readKey(plist_, domain_name, "internal energy liquid", "internal_energy_liquid");
//  dependencies_.insert(ie_key_);

};

Teuchos::RCP<Evaluator>
LakeEnthalpyEvaluator::Clone() const {
  return Teuchos::rcp(new LakeEnthalpyEvaluator(*this));
};


void LakeEnthalpyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Teuchos::OSTab tab = vo_.getOSTab();
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> u_l = S.GetPtr<CompositeVector>(ie_key_,tag);
  *result[0] = *u_l;


  if (include_work_) {
    Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_,tag);
    Teuchos::RCP<const CompositeVector> n_l = S.GetPtr<CompositeVector>(dens_key_,tag);

    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp,false);
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] += 1.e-6*pres_v[0][i]/nl_v[0][i];
      }
    }
  }
};


void LakeEnthalpyEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // not implemented
  if (wrt_key == ie_key_) {
    result[0]->PutScalar(1.);
  } else if (wrt_key ==pres_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> n_l = S.GetPtr<CompositeVector>(dens_key_,tag);

    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] = 1.e-6/nl_v[0][i];
      }
    }

  } else if (wrt_key ==dens_key_) {
    AMANZI_ASSERT(include_work_);
    
    Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_,tag);
    Teuchos::RCP<const CompositeVector> n_l = S.GetPtr<CompositeVector>(dens_key_,tag);

    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp,false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        // 1.e-6 converts to MJoules
        result_v[0][i] = -1.e-6*pres_v[0][i]/std::pow(nl_v[0][i], 2);
      }
    }
  }
};


} //namespace
} //namespace
