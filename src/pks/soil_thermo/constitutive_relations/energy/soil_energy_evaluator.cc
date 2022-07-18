/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

Evaluator for water density.
----------------------------------------------------------------------------- */


#include "soil_energy_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilEnergyEvaluator::SoilEnergyEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{temperature_key_, tag});

  // -- density
  density_key_ = Keys::readKey(plist_, domain_name, "density", "density");
  dependencies_.insert(KeyTag{density_key_, tag});

  // -- heat capacity
  heat_capacity_key_ = Keys::readKey(plist_, domain_name, "heat capacity", "heat_capacity");
  dependencies_.insert(KeyTag{heat_capacity_key_, tag});

  // -- pressure
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(KeyTag{pres_key_, tag});

};

SoilEnergyEvaluator::SoilEnergyEvaluator(const SoilEnergyEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    temperature_key_(other.temperature_key_),
    density_key_(other.density_key_),
    heat_capacity_key_(other.heat_capacity_key_),
    pres_key_(other.pres_key_){};

Teuchos::RCP<Evaluator>
SoilEnergyEvaluator::Clone() const {
  return Teuchos::rcp(new SoilEnergyEvaluator(*this));
};


void SoilEnergyEvaluator::Evaluate_(const State& S,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

  double rho0 = 1200.;
  double cp0 = 800./rho0;

  // evaluate density
  const Epetra_MultiVector& rho =
  *S.GetPtr<CompositeVector>(density_key_,tag)->ViewComponent("cell",false);

  // evaluate heat capacity
  const Epetra_MultiVector& cp =
  *S.GetPtr<CompositeVector>(heat_capacity_key_,tag)->ViewComponent("cell",false);

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i];
      result_v[0][i] = rho[0][i]*cp[0][i]*T;
//      result_v[0][i] = rho0*cp0*T;
    }
  }
};


void SoilEnergyEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;

  result[0]->PutScalar(0.);

  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

    double rho0 = 1200.;
    double cp0 = 800./rho0;

    // evaluate density
    const Epetra_MultiVector& rho =
    *S.GetPtr<CompositeVector>(density_key_,tag)->ViewComponent("cell",false);

    // evaluate heat capacity
    const Epetra_MultiVector& cp =
    *S.GetPtr<CompositeVector>(heat_capacity_key_,tag)->ViewComponent("cell",false);

    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i];
        result_v[0][i] = rho[0][i]*cp[0][i];
//        result_v[0][i] = rho0*cp0;
      }
    }
  }

  if (wrt_key == pres_key_) {
      result[0]->PutScalar(0.);
  }

};


} //namespace
} //namespace
