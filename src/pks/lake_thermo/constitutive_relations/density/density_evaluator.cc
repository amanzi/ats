/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

Evaluator for water density.
----------------------------------------------------------------------------- */


#include "density_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

DensityEvaluator::DensityEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {
      
  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;    

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temperature_key_, tag });

};

Teuchos::RCP<Evaluator>
DensityEvaluator::Clone() const {
  return Teuchos::rcp(new DensityEvaluator(*this));
};


void DensityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)   
{

  Tag tag = my_keys_.front().second;

  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

//  double rho0 = 1.;
  double rho0 = 1000.;
  double a0 = 800.969e-7;
  double a1 = 588.194e-7;
  double a2 = 811.465e-8;
  double a3 = 476.600e-10;

  double temp0 = 3.85;
  double const_ampl = 1.9549e-5;
  double const_power = 1.68;

  std::string EOS_type = plist_.get<std::string>("EOS type",
      "none");   

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    std::vector<double> rho(ncomp);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i]-273.15;
      if (T > 0.) { //water
        if (EOS_type == "UNESCO") {
          rho[i] = rho0*(1+a0+a1*T-a2*T*T+a3*T*T*T); // UNESCO
        } else if (EOS_type == "Hostetler") {
          rho[i] = rho0*(1.-const_ampl*std::pow(std::fabs(T-temp0),const_power)); // Hostetler
        } else {
          std::stringstream messagestream;
          messagestream << "Lake_Thermo PK has no density EOS: " << EOS_type;
          Errors::Message message(messagestream.str());
          Exceptions::amanzi_throw(message);
         }
      }
      else { // ice 
        rho[i] = 917.;
      }
    }

    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = rho[i];
    }

    // result_v[0][0] = rho[0];
    // for (int i=1; i!=ncomp; ++i) {
    //   result_v[0][i] = 0.5*(rho[i]+rho[i-1]);
    // }

  }
};


void DensityEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  result[0]->PutScalar(0.);

  std::string EOS_type = plist_.get<std::string>("EOS type",
      "none");   

  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

    double rho0 = 1000.;
    double a0 = 800.969e-7;
    double a1 = 588.194e-7;
    double a2 = 811.465e-8;
    double a3 = 476.600e-10;

    double temp0 = 3.85;
    double const_ampl = 1.9549e-5;
    double const_power = 1.68;

    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      std::vector<double> drho(ncomp);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i]-273.15;
        if (T > 0.) { //water
          if (EOS_type == "UNESCO") {
            drho[i] = rho0*(a1 - 2.*a2*T + 3.*a3*T*T); // UNESCO
          } else if (EOS_type == "Hostetler") {
            double x = T-temp0;
            double sign = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
            drho[i] = rho0*(-const_power*const_ampl*std::pow(std::fabs(T-temp0),const_power-1)*sign); // Hostetler
          } else {
            std::stringstream messagestream;
            messagestream << "Lake_Thermo PK has no density EOS: " << EOS_type;
            Errors::Message message(messagestream.str());
            Exceptions::amanzi_throw(message);
          }          
        }
        else { // ice 
          drho[i] = 0.;
        }
      }

      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = drho[i];
      }

      // result_v[0][0] = drho[0];
      // for (int i=1; i!=ncomp; ++i) {
      //   result_v[0][i] = 0.5*(drho[i]+drho[i-1]);
      // }

    }

  } // if
};


} //namespace
} //namespace
