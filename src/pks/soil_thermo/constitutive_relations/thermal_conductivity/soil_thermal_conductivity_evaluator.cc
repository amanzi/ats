/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity of soil model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "soil_thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilThermalConductivityEvaluator::SoilThermalConductivityEvaluator(
    Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  //  // -- temperature
  //  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  //  dependencies_.insert(temperature_key_);

  //  cell_is_ice_key_ = Keys::readKey(plist_, domain_name, "ice", "ice");
  //  dependencies_.insert(cell_is_ice_key_);

  //  // -- cell volume
  //  cell_vol_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  //  dependencies_.insert(cell_vol_key_);

  // -- water content
  water_content_key_ = Keys::readKey(plist_, domain_name, "water content", "water_content");
  dependencies_.insert(KeyTag{water_content_key_, tag});

  // -- ice content
  ice_content_key_ = Keys::readKey(plist_, domain_name, "ice content", "ice_content");
  dependencies_.insert(KeyTag{ice_content_key_, tag});

  // -- ice content
  pressure_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(KeyTag{pressure_key_, tag});

  //  AMANZI_ASSERT(plist_.isSublist("soil thermal conductivity parameters"));
  //  Teuchos::ParameterList sublist = plist_.sublist("soil thermal conductivity parameters");

}


SoilThermalConductivityEvaluator::SoilThermalConductivityEvaluator(
    const SoilThermalConductivityEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
        temperature_key_(other.temperature_key_),
        water_content_key_(other.water_content_key_),
        ice_content_key_(other.ice_content_key_),
        pressure_key_(other.pressure_key_){}


Teuchos::RCP<Evaluator>
SoilThermalConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new SoilThermalConductivityEvaluator(*this));
}

void SoilThermalConductivityEvaluator::Evaluate_(const State& S,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;

  //  // get temperature
  //  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

  // get water content
  Teuchos::RCP<const CompositeVector> wc = S.GetPtr<CompositeVector>(water_content_key_,tag);

  //  // get ice content
  //  Teuchos::RCP<const CompositeVector> ic = S.GetPtr<CompositeVector>(ice_content_key_,tag);

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();

  double eps = 1.e-10;

  for (CompositeVector::name_iterator comp=result[0]->begin();
      comp!=result[0]->end(); ++comp) {
    // much more efficient to pull out vectors first
    //      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp,false);
    //      const Epetra_MultiVector& ic_v = *ic->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);

    for (int i=0; i!=ncomp; ++i) {

      // read these parameters from states
      double por = 0.;
      //        double wl = 0.;
      double wi = 0.;
      double wl = wc_v[0][i];// * 1.8e-5; // CONVERTED UNITS
      //        double wi = ic_v[0][i];

      //        std::cout << "soil water content wc_v[0][i] = " << wc_v[0][i] << std::endl;
      //        std::cout << "wl = " << wl << std::endl;

      double lambda_sat;
      double lambda_dry;

      double row0, roi, row0_d_roi, roi_d_row0, lamw0, lami;

      double rosdry = 1200.;

      double CK_const, CK_const1, CK_const2;

      row0  = 1000.;    // reference density of water,    kg/m**3
      roi   = 917.;     // density of ice,    kg/m**3

      lami  = 2.2;      // molecular conductivity of ice, J/(m*s*K)
      lamw0 = 0.561;    // molecular conductivity of water,     J/(m*s*K)

      row0_d_roi = row0/roi;
      roi_d_row0 = roi/row0;

      double CK_consts[4], CK_consts1[3], CK_consts2[3];

      CK_consts[0] = 4.6;  // for gravel and coarse sand
      CK_consts[1] = 3.55; // for medium and fine sand
      CK_consts[2] = 1.9;  // silty and clay soils
      CK_consts[3] = 0.6;  // organic fibrous soils

      CK_consts1[0] = 1.7;  // for crashed rocks
      CK_consts1[1] = 0.75; // for mineral soils
      CK_consts1[2] = 0.3;  // organic fibrous soils

      CK_consts2[0] = 1.8;  // for crashed rocks
      CK_consts2[1] = 1.2;  // for mineral soils
      CK_consts2[2] = 0.87; // organic fibrous soils

      double quartz_ratio = 0.1;
      double lambda_quartz = 7.7;
      double lambda_othmin;
      if (quartz_ratio > 0.2) {
        lambda_othmin = 2.;
      }
      else {
        lambda_othmin = 3.;
      }
      double lambda_solids = pow(lambda_quartz,quartz_ratio) * pow(lambda_othmin,(1.-quartz_ratio));

      // Conversion from mass ratios to volume ratios
      double water_vol_ratio = wl / (por*(wl + row0/roi*wi + row0/rosdry) + eps);
      double ice_vol_ratio = wi / (por*(wi + roi/row0*wl + roi/rosdry) + eps);

      double waterice_vol_ratio = water_vol_ratio + ice_vol_ratio;

      CK_const = CK_consts[2]; // silty and clay soils are assumed
      double Kersten = CK_const*waterice_vol_ratio / (1. + (CK_const - 1.)*waterice_vol_ratio + eps);

      double water_sat_ratio = por*water_vol_ratio/(waterice_vol_ratio + eps);
      double ice_sat_ratio = por*ice_vol_ratio/(waterice_vol_ratio + eps);

      // This is the original formula from Johansen (1975) extended for the case
      // with the ice content
      lambda_sat = pow(lambda_solids,(1. - por)) * pow(lamw0,(water_sat_ratio))*pow(lami,(ice_sat_ratio));

      CK_const1 = CK_consts1[1]; // mineral soils are assumed
      CK_const2 = CK_consts2[1]; // mineral soils are assumed
      // Cote and Konrad (2005) formula for heat conduction coefficient of dry soil
      lambda_dry = CK_const1*pow(10.,(-CK_const2*por));

      //        result_v[0][i] = 0.75;
      result_v[0][i] = (lambda_sat - lambda_dry)*Kersten + lambda_dry;


    } // i
  }

}


void SoilThermalConductivityEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& result) {
  Tag tag = my_keys_.front().second;
  //  std::cout<<"SOIL THERMAL CONDUCITIVITY: Derivative not implemented yet!"<<wrt_key<<"\n";
  //  AMANZI_ASSERT(0); // not implemented, not yet needed
  //  result[0]->Scale(1.e-6); // convert to MJ
  result[0]->PutScalar(0.);
}

} //namespace
} //namespace
