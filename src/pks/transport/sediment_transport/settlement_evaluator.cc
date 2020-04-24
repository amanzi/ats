/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors:
*/

#include "settlement_evaluator.hh"
#include "boost/math/constants/constants.hpp"

namespace Amanzi {

SettlementRateEvaluator :: SettlementRateEvaluator(Teuchos::ParameterList& plist) :
  SecondaryVariableFieldEvaluator(plist) {

  Key domain_name = "surface";
  
  velocity_key_ = plist_.get<std::string>("velocity key",
                                     Keys::getKey(domain_name,"velocity"));

  sediment_key_ =  plist_.get<std::string>("sediment key",
                                     Keys::getKey(domain_name,"sediment"));

  tau_d_ = plist_.get<double>("critical shear stress");
  ws_ = plist_.get<double>("settling velocity");
  gamma_ = plist_.get<double>("specific weight of water");
  sediment_density_ = plist_.get<double>("sediment density [kg/m^3]");

  
  umax_ = plist_.get<double>("max current");
  xi_ = plist_.get<double>("Chezy parameter");

  double pi = boost::math::constants::pi<double>();

  lambda_ = 8./(3*pi) * (umax_/(xi_*xi_));
    
  dependencies_.insert("surface-pressure");
  dependencies_.insert(sediment_key_);
    
}

  
SettlementRateEvaluator ::SettlementRateEvaluator (const SettlementRateEvaluator & other) :
  SecondaryVariableFieldEvaluator(other),
  velocity_key_(other.velocity_key_), sediment_key_(other.sediment_key_) {

  tau_d_ = other.tau_d_;
  ws_ = other.ws_;
  gamma_ = other.gamma_;
  lambda_ = other.lambda_;
  sediment_density_ = other.sediment_density_;
} 


Teuchos::RCP<FieldEvaluator> SettlementRateEvaluator ::Clone() const {
  return Teuchos::rcp(new SettlementRateEvaluator (*this));
}


void SettlementRateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  const Epetra_MultiVector& vel = *S->GetFieldData(velocity_key_)->ViewComponent("cell");
  const Epetra_MultiVector& tcc = *S->GetFieldData(sediment_key_)->ViewComponent("cell");
  Epetra_MultiVector& result_c = *result->ViewComponent("cell");
  
  for (int c=0; c<result_c.MyLength(); c++){
    double tau_0 = gamma_ * lambda_ * (sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]));
    
    if (tau_0 < tau_d_){
      result_c[0][c] = sediment_density_ * ws_ * std::min(tcc[0][c], 0.5) * (1 - tau_0 / tau_d_);
    }else{
      result_c[0][c] = 0.;
    }

  }
   
}

void SettlementRateEvaluator::EvaluateFieldPartialDerivative_ (const Teuchos::Ptr<State>& S,
                                                            Key wrt_key,
                                                            const Teuchos::Ptr<CompositeVector>& result) {
   AMANZI_ASSERT(0); 
}
  
  
} // namespace
