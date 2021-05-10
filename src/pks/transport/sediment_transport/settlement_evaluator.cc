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
  EvaluatorSecondaryMonotypeCV(plist) {

  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  
  velocity_key_ = plist_.get<std::string>("velocity key",
                                     Keys::getKey(domain_name,"velocity"));

  sediment_key_ =  plist_.get<std::string>("sediment key",
                                     Keys::getKey(domain_name,"sediment"));

  tau_d_ = plist_.get<double>("critical shear stress");
  ws_ = plist_.get<double>("settling velocity");
  gamma_ = plist_.get<double>("specific weight of water");
  sediment_density_ = plist_.get<double>("sediment density [kg m^-3]");

  
  umax_ = plist_.get<double>("max current");
  xi_ = plist_.get<double>("Chezy parameter");
  Cf_ = plist_.get<double>("drag coefficient");

  double pi = boost::math::constants::pi<double>();

  lambda_ = 8./(3*pi) * (umax_/(xi_*xi_));
    
  dependencies_.insert(KeyTag{"surface-pressure", tag});
  dependencies_.insert(KeyTag{sediment_key_, tag});
    
}

  
SettlementRateEvaluator ::SettlementRateEvaluator (const SettlementRateEvaluator & other) :
  EvaluatorSecondaryMonotypeCV(other),
  velocity_key_(other.velocity_key_), sediment_key_(other.sediment_key_) {

  tau_d_ = other.tau_d_;
  ws_ = other.ws_;
  gamma_ = other.gamma_;
  lambda_ = other.lambda_;
  sediment_density_ = other.sediment_density_;
  Cf_ = other.Cf_;
} 


Teuchos::RCP<Evaluator> SettlementRateEvaluator ::Clone() const {
  return Teuchos::rcp(new SettlementRateEvaluator (*this));
}


void SettlementRateEvaluator::Evaluate_(const State& S,
                                        const std::vector<CompositeVector*>& result){

  Tag tag = my_keys_.front().second;
  
  const Epetra_MultiVector& vel = *S.GetPtr<CompositeVector>(velocity_key_, tag)->ViewComponent("cell");
  const Epetra_MultiVector& tcc = *S.GetPtr<CompositeVector>(sediment_key_, tag)->ViewComponent("cell");
  Epetra_MultiVector& result_c = *result[0]->ViewComponent("cell");
  
  for (int c=0; c<result_c.MyLength(); c++){
    //double tau_0 = gamma_ * lambda_ * (sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]));
    double tau_0 = gamma_ * Cf_ * (sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c])*sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]));
    
    if (tau_0 < tau_d_){
      result_c[0][c] = sediment_density_ * ws_ * std::min(tcc[0][c], 0.5) * (1 - tau_0 / tau_d_);
    }else{
      result_c[0][c] = 0.;
    }

  }
   
}

void SettlementRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key, const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& result){
   AMANZI_ASSERT(0); 
}
  
  
} // namespace
