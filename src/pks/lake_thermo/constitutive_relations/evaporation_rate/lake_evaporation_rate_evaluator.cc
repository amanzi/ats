/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a heat capacity of soil model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "lake_evaporation_rate_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

LakeEvaporationRateEvaluator::LakeEvaporationRateEvaluator(
    Teuchos::ParameterList& plist) :
            SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("lake evaporation rate key",
        "surface-evaporation_rate");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

}

LakeEvaporationRateEvaluator::LakeEvaporationRateEvaluator(
    const LakeEvaporationRateEvaluator& other) :
            SecondaryVariableFieldEvaluator(other),
            temperature_key_(other.temperature_key_){}


Teuchos::RCP<FieldEvaluator>
LakeEvaporationRateEvaluator::Clone() const {
  return Teuchos::rcp(new LakeEvaporationRateEvaluator(*this));
}

void LakeEvaporationRateEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) {

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  // read parameters from the met data
  Teuchos::ParameterList& param_list = plist_.sublist("parameters");
  Amanzi::FunctionFactory fac;
  Teuchos::RCP<Amanzi::Function> SS_func_ = Teuchos::rcp(fac.Create(param_list.sublist("solar radiation")));
  Teuchos::RCP<Amanzi::Function> E_a_func_ = Teuchos::rcp(fac.Create(param_list.sublist("atmospheric downward radiation")));
  Teuchos::RCP<Amanzi::Function> T_a_func_ = Teuchos::rcp(fac.Create(param_list.sublist("air temperature")));
  Teuchos::RCP<Amanzi::Function> H_a_func_ = Teuchos::rcp(fac.Create(param_list.sublist("air humidity")));
  Teuchos::RCP<Amanzi::Function> P_a_func_ = Teuchos::rcp(fac.Create(param_list.sublist("atmospheric pressure")));

  std::vector<double> args(1);
  args[0] = S->time();
  double SS = (*SS_func_)(args);
  double E_a = (*E_a_func_)(args);
  double T_a = (*T_a_func_)(args);
  double q_a = (*H_a_func_)(args);
  double P_a = (*P_a_func_)(args);

  for (CompositeVector::name_iterator comp=result->begin();
      comp!=result->end(); ++comp) {

    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);

    // get temperature
    const Epetra_MultiVector& temp_v = *S->GetFieldData(temperature_key_)
          ->ViewComponent("cell",false);

    double T_s = temp_v[0][ncomp-1];

    double b1_vap   = 610.78;        // Coefficient [N m^{-2} = kg m^{-1} s^{-2}]
    double b3_vap   = 273.16;        // Triple point [K]
    double b2w_vap  = 17.2693882;    // Coefficient (water)
    double b2i_vap  = 21.8745584;    // Coefficient (ice)
    double b4w_vap  = 35.86;         // Coefficient (temperature) [K]
    double b4i_vap  = 7.66;          // Coefficient (temperature) [K]

    // Saturation water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
    double wvpres_s;

    double h_ice = 0.;
    double h_Ice_min_flk = 1.e-9;
    if (h_ice < h_Ice_min_flk) {  // Water surface
      wvpres_s = b1_vap*std::exp(b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap));
    } else {                      // Ice surface
      wvpres_s = b1_vap*std::exp(b2i_vap*(T_s-b3_vap)/(T_s-b4i_vap));
    }

    // Saturation specific humidity at T=T_s
    double tpsf_R_dryair    = 2.8705e2;  // Gas constant for dry air [J kg^{-1} K^{-1}]
    double tpsf_R_watvap    = 4.6151e2;  // Gas constant for water vapour [J kg^{-1} K^{-1}]
    double tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap;
    double q_s = tpsf_Rd_o_Rv*wvpres_s/(P_a-(1.-tpsf_Rd_o_Rv)*wvpres_s);

    double height_tq = 2.;
    double tpsf_kappa_q_a   = 2.4e-05; // Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

    double LE = -tpsf_kappa_q_a*(q_a-q_s)/height_tq;

    double rho_a = P_a/tpsf_R_dryair/T_s/(1.+(1./tpsf_Rd_o_Rv-1.)*q_s);

    double tpsf_c_a_p  = 1.005e3; // Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
    double tpsf_L_evap = 2.501e6; // Specific heat of evaporation [J kg^{-1}]
    double tpl_L_f     = 3.3e5;   // Latent heat of fusion [J kg^{-1}]

    double Q_watvap   = LE*rho_a;
    LE = tpsf_L_evap;
    if (h_ice >= h_Ice_min_flk) LE = LE + tpl_L_f;   // Add latent heat of fusion over ice
    LE = Q_watvap*LE;

    double row0 = 1.e+3;
    double evap_rate = LE/(row0*tpsf_L_evap); //*row0;
    evap_rate = abs(evap_rate);

    for (int i=0; i!=ncomp; ++i) {

      result_v[0][i] = evap_rate;

    } // i
  }

}


void LakeEvaporationRateEvaluator::EvaluateFieldPartialDerivative_(
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
