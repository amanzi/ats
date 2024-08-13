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
    EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  std::cout << "surf temp eval domain_name = " << domain_name << std::endl;

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temperature_key_, tag });

  std::cout << "temperature_key_ = " << temperature_key_ << std::endl;

};

Teuchos::RCP<Evaluator>
LakeSurfaceTemperatureEvaluator::Clone() const {
  return Teuchos::rcp(new LakeSurfaceTemperatureEvaluator(*this));
};

void LakeSurfaceTemperatureEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  std::cout << "In surface temp eval EnsureCompatibility" << std::endl;

}


void LakeSurfaceTemperatureEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  std::cout << "SurfTempEval check 1 " << std::endl;
  std::cout << "temperature_key_ = " << temperature_key_ << std::endl;

  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);
  int ncomp_temp = temp->size("cell", false);
  const Epetra_MultiVector& temp_v = *temp->ViewComponent("cell",false);

  double T_surf = temp_v[0][ncomp_temp-1];

  std::cout << "T_surf = " << T_surf << std::endl;

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = T_surf;
    }
  }
};


void LakeSurfaceTemperatureEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  result[0]->PutScalar(0.);

};


} //namespace
} //namespace
