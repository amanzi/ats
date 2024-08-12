/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a heat flux at the surface of lake model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "heat_flux_bc_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

HeatFluxBCEvaluator::HeatFluxBCEvaluator(
    Teuchos::ParameterList& plist) :
            SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("heat flux bc key",
        "surface-heat_flux_bc");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  // -- thermal conductivity
  conductivity_key_ = Keys::readKey(plist_, domain_name, "thermal conductivity", "thermal_conductivity");
  dependencies_.insert(conductivity_key_);

  //  AMANZI_ASSERT(plist_.isSublist("heat flux bc parameters"));
  //  Teuchos::ParameterList sublist = plist_.sublist("heat flux bc parameters");

  // later: read these parameters from xml
  SS = 0.;      // solar radiation (read from met data)
  alpha_w = 0.06; // water albedo
  alpha_i = 0.40; // ice albedo
  E_a = 0.;     // atmospheric downward radiation (read from met data)
  E_s = 0.;     // surface radiation (Stefan-Boltzman law)
  H = 0.;       // "sensible" heat
  LE = 0.;      // latent heat

}


HeatFluxBCEvaluator::HeatFluxBCEvaluator(
    const HeatFluxBCEvaluator& other) :
            SecondaryVariableFieldEvaluator(other),
            SS(other.SS),
            alpha_w(other.alpha_w),
            alpha_i(other.alpha_i),
            E_a(other.E_a),
            E_s(other.E_s),
            H(other.H),
            LE(other.LE),
            temperature_key_(other.temperature_key_),
            conductivity_key_(other.temperature_key_){}


Teuchos::RCP<FieldEvaluator>
HeatFluxBCEvaluator::Clone() const {
  return Teuchos::rcp(new HeatFluxBCEvaluator(*this));
}

void HeatFluxBCEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) {

  ice_cover_ = false; // first always assume that there is no ice

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  // get conductivity
  Teuchos::RCP<const CompositeVector> cond = S->GetFieldData(conductivity_key_);

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
  SS = (*SS_func_)(args);
  E_a = (*E_a_func_)(args);
  double T_a = (*T_a_func_)(args);
  double q_a = (*H_a_func_)(args);
  double P_a = (*P_a_func_)(args);

  double sigma = 5.67e-8;     // Stefan-Boltzman constant
  double c_lwrad_emis = 0.99; //Surface emissivity with respect to the long-wave radiation

  for (CompositeVector::name_iterator comp=result->begin();
      comp!=result->end(); ++comp) {

    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    const Epetra_MultiVector& cond_v = *cond->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);

    for (int i=0; i!=ncomp; ++i) {

      double T_s = temp_v[0][ncomp-1];

      E_s = c_lwrad_emis*sigma*pow(T_s,4);

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
      double tpsf_kappa_t_a   = 2.2e-05; // Molecular temperature conductivity of air [m^{2} s^{-1}]
      double tpsf_kappa_q_a   = 2.4e-05; // Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

      H = -tpsf_kappa_t_a*(T_a-T_s)/height_tq;
      LE = -tpsf_kappa_q_a*(q_a-q_s)/height_tq;

      double rho_a = P_a/tpsf_R_dryair/T_s/(1.+(1./tpsf_Rd_o_Rv-1.)*q_s);

      double tpsf_c_a_p  = 1.005e3; // Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
      double tpsf_L_evap = 2.501e6; // Specific heat of evaporation [J kg^{-1}]
      double tpl_L_f     = 3.3e5;   // Latent heat of fusion [J kg^{-1}]

      H = H*rho_a*tpsf_c_a_p;
      double Q_watvap   = LE*rho_a;
      LE = tpsf_L_evap;
      if (h_ice >= h_Ice_min_flk) LE = LE + tpl_L_f;   // Add latent heat of fusion over ice
      LE = Q_watvap*LE;

     double alpha;
     alpha = (T_s < 273.15) ? alpha_i : alpha_w;

     double hour_sec = 60.*60;
     double interval = 24.;

     int n = int(S->time());
     int day = n / (24 * 3600);
     n = n % (24 * 3600);
     int hour = n / 3600;
     n %= 3600;
     int minutes = n / 60 ;
     n %= 60;
     int seconds = n;

    //  std::cout << day << " " << "days " << hour
    //      << " " << "hours " << minutes << " "
    //      << "minutes " << seconds << " "
    //      << "seconds "  << std::endl;

     double pi = 3.1415;
     double longitude = -164.456696;
     double latitude = 66.558765;
     double phi   = latitude*pi/180.0;
     double delta = longitude*pi/180.0; //23.5*pi/180.0*cos(2*pi*(double(day)-173.0)/365.0);
     double theta = pi*(hour-12.0)/12.0;
     double sinh0 = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(theta);
     sinh0 = fmax(sinh0,0.0); 

    //  std::cout << "sinh0 = " << sinh0 << std::endl;

    //  SS *= sqrt(1.-sinh0*sinh0);

    //  SS = SS/(hour_sec*interval);
      // SS = 0.92*SS; // FoxDen
      SS = 0.9*SS; // Atqasuk
      // SS = 0.9*SS; // Toolik

      // SS = (T_s < 273.15) ? 0. : SS;

      // std::cout << "SS = " << SS << ", E_a = " << E_a << ", E_s = " << E_s << ", H = " << H << ", LE = " << LE << std::endl;
      
      result_v[0][i] = SS*(1.-alpha) + E_a - E_s - H - LE; 

      result_v[0][i] *= -1.; ///cond_v[0][i];

      if (isnan(result_v[0][i])) {
          std::cout << "SS = " << SS << ", E_a = " << E_a << ", E_s = " << E_s << ", H = " << H << ", LE = " << LE << std::endl;
          std::cout << "rho_a = " << rho_a << std::endl;
          std::cout << "T_s = " << T_s << std::endl;
          std::cout << "q_s = " << q_s << std::endl;
          std::cout << "tpsf_Rd_o_Rv = " << tpsf_Rd_o_Rv << std::endl;
          std::cout << "wvpres_s = " << wvpres_s << std::endl;
          std::cout << "P_a = " << P_a << std::endl;
          std::cout << "tpsf_Rd_o_Rv = " << tpsf_Rd_o_Rv << std::endl;
          std::cout << b1_vap << " " << std::exp(b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap)) << std::endl;
          std::cout << b2w_vap << std::endl;
          std::cout << T_s-b3_vap << std::endl;
          std::cout << T_s-b4w_vap << std::endl;
          std::cout << (T_s-b3_vap)/(T_s-b4w_vap) << std::endl;
          std::cout << "b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap) = " << b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap) << std::endl;
          std::cout << b1_vap << " " << std::exp(b2i_vap*(T_s-b3_vap)/(T_s-b4i_vap)) << std::endl;
          exit(0);
      }

    } // i

  }

}

void HeatFluxBCEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"HEAT FLUX BC: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
