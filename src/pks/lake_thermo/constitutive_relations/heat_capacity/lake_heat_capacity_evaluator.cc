/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a heat capacity of soil model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "lake_heat_capacity_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

LakeHeatCapacityEvaluator::LakeHeatCapacityEvaluator(
    Teuchos::ParameterList& plist) :
        EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temperature_key_, tag });

//  // -- water content
//  water_content_key_ = Keys::readKey(plist_, domain_name, "water content", "water_content");
// dependencies_.insert(KeyTag{ water_content_key_, tag });

  //  // -- ice content
  //  ice_content_key_ = Keys::readKey(plist_, domain_name, "soil ice content", "soil_ice_content");
  // dependencies_.insert(KeyTag{ ice_content_key_, tag });

  //  AMANZI_ASSERT(plist_.isSublist("soil heat capacity parameters"));
  //  Teuchos::ParameterList sublist = plist_.sublist("soil heat capacity parameters");

  double row  = 1000.; // density of water
  double roi  = 917.;  // density of ice

  // cw    = 3990.; ///row;    // specific heat of water
  cw    = 4182.; ///row; // [J/kg K]
  ci    = 2150.; ///roi;    // specific heat of ice
  // ci = cw;
}

Teuchos::RCP<Evaluator>
LakeHeatCapacityEvaluator::Clone() const {
  return Teuchos::rcp(new LakeHeatCapacityEvaluator(*this));
}

void
LakeHeatCapacityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // get temperature
  Teuchos::RCP<const CompositeVector> T = S.GetPtr<CompositeVector>(temperature_key_,tag);

//  // get water content
//  Teuchos::RCP<const CompositeVector> wc = S.GetPtr<CompositeVector>(water_content_key_,tag);

  //  // get ice content
  //  Teuchos::RCP<const CompositeVector> ic = S.GetPtr<CompositeVector>(ice_content_key_,tag);

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();

  for (CompositeVector::name_iterator comp=result[0]->begin();
      comp!=result[0]->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& T_v = *T->ViewComponent(*comp,false);
//    const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp,false);
//    const Epetra_MultiVector& ic_v = *ic->ViewComponent(*comp,false);

    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);

    std::vector<double> cc(ncomp);

    for (int i=0; i!=ncomp; ++i) {

      double T = T_v[0][i];
      cc[i] = (T < 273.15) ? ci : cw;
      // if (T < 273.15) {
      //   cc[i] = 0.5*(ci + cc[i]);
      // }
      // else {
      //   cc[i] = 0.5*(cw + cc[i]);
      // }

    } // i

    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = cc[i];
    }

    // result_v[0][0] = cc[0];
    // for (int i=1; i!=ncomp; ++i) {
    //   result_v[0][i] = 0.5*(cc[i]+cc[i-1]);
    // }

    // result_v[0][0] = cc[0];
    // result_v[0][ncomp-1] = cc[ncomp-1];
    // for (int i=1; i!=ncomp-1; ++i) {
    //   result_v[0][i] = 1./3.*(cc[i-1]+cc[i]+cc[i+1]);
    // }


  }

}

void
LakeHeatCapacityEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  result[0]->PutScalar(0.0);

//  if (wrt_key == water_content_key_) {
//
//    for (CompositeVector::name_iterator comp=result[0]->begin();
//        comp!=result[0]->end(); ++comp) {
//      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);
//
//      int ncomp = result[0]->size(*comp, false);
//      for (int i=0; i!=ncomp; ++i) {
//        result_v[0][i] = cw * 1.8e-5;
//      }
//    }
//  }
//  if (wrt_key == ice_content_key_) {
//
//    for (CompositeVector::name_iterator comp=result[0]->begin();
//        comp!=result[0]->end(); ++comp) {
//      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);
//
//      int ncomp = result[0]->size(*comp, false);
//      for (int i=0; i!=ncomp; ++i) {
//        result_v[0][i] = ci * 1.8e-5;
//      }
//    }
//  }

}

} //namespace
} //namespace
