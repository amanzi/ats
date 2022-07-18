/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a heat flux at the surface of soil model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "soil_heat_flux_bc_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilHeatFluxBCEvaluator::SoilHeatFluxBCEvaluator(
    Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{temperature_key_, tag});

  AMANZI_ASSERT(plist_.isSublist("soil heat flux bc parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("soil heat flux bc parameters");

  // later: read these parameters from xml
  SS = 0;
  alpha = 0;
  E_a = 0;
  E_s = 0;
  H = 0;
  LE = 0;

}


SoilHeatFluxBCEvaluator::SoilHeatFluxBCEvaluator(
    const SoilHeatFluxBCEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
        SS(other.SS),
        alpha(other.alpha),
        E_a(other.E_a),
        E_s(other.E_s),
        H(other.H),
        LE(other.LE),
        temperature_key_(other.temperature_key_){}


Teuchos::RCP<Evaluator>
SoilHeatFluxBCEvaluator::Clone() const {
  return Teuchos::rcp(new SoilHeatFluxBCEvaluator(*this));
}

void SoilHeatFluxBCEvaluator::Evaluate_(const State& S,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;

  ice_cover_ = false; // first always assume that there is no ice

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

  for (CompositeVector::name_iterator comp=result[0]->begin();
      comp!=result[0]->end(); ++comp) {

    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);

    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = SS*(1.-alpha) + E_a - E_s - H - LE;
    } // i

  }

}

void SoilHeatFluxBCEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;
  std::cout<<"HEAT FLUX BC: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result[0]->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
