/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a heat capacity of soil model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "soil_heat_capacity_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilHeatCapacityEvaluator::SoilHeatCapacityEvaluator(
    Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  //  // -- temperature
  //  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  //  dependencies_.insert(temperature_key_);

  // -- water content
  water_content_key_ = Keys::readKey(plist_, domain_name, "water content", "water_content");
  dependencies_.insert(KeyTag{water_content_key_, tag});

  //  // -- ice content
  //  ice_content_key_ = Keys::readKey(plist_, domain_name, "soil ice content", "soil_ice_content");
  //  dependencies_.insert(ice_content_key_);

  //  AMANZI_ASSERT(plist_.isSublist("soil heat capacity parameters"));
  //  Teuchos::ParameterList sublist = plist_.sublist("soil heat capacity parameters");

  double row  = 1000.; // density of water
  double roi  = 917.;  // density of ice
  double ros  = 1200.; // density of soil

  cw    = 3990./row;    // specific heat of water
  ci    = 2150./roi;    // specific heat of ice
  cg    = 800.; ///ros;     // specific heat of dry soils

}


SoilHeatCapacityEvaluator::SoilHeatCapacityEvaluator(
    const SoilHeatCapacityEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
        cg(other.cg),
        cw(other.cw),
        ci(other.ci),
        water_content_key_(other.water_content_key_),
        ice_content_key_(other.ice_content_key_){}


Teuchos::RCP<Evaluator>
SoilHeatCapacityEvaluator::Clone() const
{
  return Teuchos::rcp(new SoilHeatCapacityEvaluator(*this));
}

void SoilHeatCapacityEvaluator::Evaluate_(const State& S,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;

  // get water content
  Teuchos::RCP<const CompositeVector> wc = S.GetPtr<CompositeVector>(water_content_key_,tag);

  //  // get ice content
  //  Teuchos::RCP<const CompositeVector> ic = S.GetPtr<CompositeVector>(ice_content_key_,tag);

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();

  for (CompositeVector::name_iterator comp=result[0]->begin();
      comp!=result[0]->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp,false);
    //      const Epetra_MultiVector& ic_v = *ic->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);

    for (int i=0; i!=ncomp; ++i) {

      double W = wc_v[0][i]; // * 1.8e-5; // CONVERTED UNITS from mol/m^3 to volume ratio
      double I = 0.;//= ic_v[0][i] * 1.8e-5;

      result_v[0][i] = cg + cw*W + ci*I;

      //          std::cout << "W = " << W << std::endl;
      //          std::cout << "i = " << i << ", heat capacity = " << result_v[0][i] << std::endl;

    } // i
  }

}


void SoilHeatCapacityEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;
  result[0]->PutScalar(0.0);
  if (wrt_key == water_content_key_) {

    for (CompositeVector::name_iterator comp=result[0]->begin();
        comp!=result[0]->end(); ++comp) {
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = cw * 1.8e-5;
      }
    }
  }
  if (wrt_key == ice_content_key_) {

    for (CompositeVector::name_iterator comp=result[0]->begin();
        comp!=result[0]->end(); ++comp) {
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = ci * 1.8e-5;
      }
    }
  }

}

} //namespace
} //namespace
