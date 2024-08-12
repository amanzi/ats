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
        SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("lake heat capacity key",
        "surface-heat_capacity");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

//  // -- water content
//  water_content_key_ = Keys::readKey(plist_, domain_name, "water content", "water_content");
//  dependencies_.insert(water_content_key_);

  //  // -- ice content
  //  ice_content_key_ = Keys::readKey(plist_, domain_name, "soil ice content", "soil_ice_content");
  //  dependencies_.insert(ice_content_key_);

  //  AMANZI_ASSERT(plist_.isSublist("soil heat capacity parameters"));
  //  Teuchos::ParameterList sublist = plist_.sublist("soil heat capacity parameters");

  double row  = 1000.; // density of water
  double roi  = 917.;  // density of ice

  // cw    = 3990.; ///row;    // specific heat of water
  cw    = 4182.; ///row; // [J/kg K]
  ci    = 2150.; ///roi;    // specific heat of ice
  // ci = cw;
}


LakeHeatCapacityEvaluator::LakeHeatCapacityEvaluator(
    const LakeHeatCapacityEvaluator& other) :
        SecondaryVariableFieldEvaluator(other),
        cw(other.cw),
        ci(other.ci),
        temperature_key_(other.temperature_key_),
        water_content_key_(other.water_content_key_),
        ice_content_key_(other.ice_content_key_){}


Teuchos::RCP<FieldEvaluator>
LakeHeatCapacityEvaluator::Clone() const {
  return Teuchos::rcp(new LakeHeatCapacityEvaluator(*this));
}

void LakeHeatCapacityEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) {

  // get temperature
  Teuchos::RCP<const CompositeVector> T = S->GetFieldData(temperature_key_);

//  // get water content
//  Teuchos::RCP<const CompositeVector> wc = S->GetFieldData(water_content_key_);

  //  // get ice content
  //  Teuchos::RCP<const CompositeVector> ic = S->GetFieldData(ice_content_key_);

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  for (CompositeVector::name_iterator comp=result->begin();
      comp!=result->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& T_v = *T->ViewComponent(*comp,false);
//    const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp,false);
//    const Epetra_MultiVector& ic_v = *ic->ViewComponent(*comp,false);

    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);

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


void LakeHeatCapacityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  result->PutScalar(0.0);

//  if (wrt_key == water_content_key_) {
//
//    for (CompositeVector::name_iterator comp=result->begin();
//        comp!=result->end(); ++comp) {
//      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);
//
//      int ncomp = result->size(*comp, false);
//      for (int i=0; i!=ncomp; ++i) {
//        result_v[0][i] = cw * 1.8e-5;
//      }
//    }
//  }
//  if (wrt_key == ice_content_key_) {
//
//    for (CompositeVector::name_iterator comp=result->begin();
//        comp!=result->end(); ++comp) {
//      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);
//
//      int ncomp = result->size(*comp, false);
//      for (int i=0; i!=ncomp; ++i) {
//        result_v[0][i] = ci * 1.8e-5;
//      }
//    }
//  }

}

} //namespace
} //namespace
