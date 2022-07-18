/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

Evaluator for soil density.
----------------------------------------------------------------------------- */


#include "soil_density_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilDensityEvaluator::SoilDensityEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{temperature_key_, tag});

};

SoilDensityEvaluator::SoilDensityEvaluator(const SoilDensityEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    temperature_key_(other.temperature_key_) {};

Teuchos::RCP<Evaluator>
SoilDensityEvaluator::Clone() const {
  return Teuchos::rcp(new SoilDensityEvaluator(*this));
}


void SoilDensityEvaluator::Evaluate_(const State& S,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;

  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

  double rho0 = 1200.;

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i];
      result_v[0][i] = rho0; // * (1.+8.0*1.e-5 + 5.88*1.e-5*T - 8.11*1.e-6*T*T + 4.77*1.e-8*T*T*T);
    }
  }
};


void SoilDensityEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;
  // not implemented
  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

    double rho0 = 1200.;

    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i];
        result_v[0][i] = 0.; //rho0 * (5.88*1.e-5 - 2.*8.11*1.e-6*T + 3.*4.77*1.e-8*T*T);
      }
    }
  }
};


} //namespace
} //namespace
