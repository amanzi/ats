/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for surface temperature.
----------------------------------------------------------------------------- */


#include "lake_surface_temperature_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

LakeSurfaceTemperatureEvaluator::LakeSurfaceTemperatureEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_.empty()) {

    my_key_ = plist_.get<std::string>("temperature", "temperature");
  }

  std::cout << "my_key_ = " << my_key_ << std::endl;

  // Set up my dependencies.
  std::string domain_name = "domain"; //Keys::getDomain(my_key_);

  std::cout << "surf temp eval domain_name = " << domain_name << std::endl;

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  std::cout << "temperature_key_ = " << temperature_key_ << std::endl;

};

LakeSurfaceTemperatureEvaluator::LakeSurfaceTemperatureEvaluator(const LakeSurfaceTemperatureEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    temperature_key_(other.temperature_key_){};

Teuchos::RCP<FieldEvaluator>
LakeSurfaceTemperatureEvaluator::Clone() const {
  return Teuchos::rcp(new LakeSurfaceTemperatureEvaluator(*this));
};

void LakeSurfaceTemperatureEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  std::cout << "In surface temp eval EnsureCompatibility" << std::endl;

}


void LakeSurfaceTemperatureEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  std::cout << "SurfTempEval check 1 " << std::endl;
  std::cout << "temperature_key_ = " << temperature_key_ << std::endl;

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);
  int ncomp_temp = temp->size("cell", false);
  const Epetra_MultiVector& temp_v = *temp->ViewComponent("cell",false);

  double T_surf = temp_v[0][ncomp_temp-1];

  std::cout << "T_surf = " << T_surf << std::endl;

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = T_surf;
    }
  }
};


void LakeSurfaceTemperatureEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  result->PutScalar(0.);

};


} //namespace
} //namespace
