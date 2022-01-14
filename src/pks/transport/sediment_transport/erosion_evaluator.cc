/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "erosion_evaluator.hh"
#include "boost/math/constants/constants.hpp"

namespace Amanzi {

ErosionRateEvaluator :: ErosionRateEvaluator(Teuchos::ParameterList& plist) :
  EvaluatorSecondaryMonotypeCV(plist) {

  Key domain_name = "surface";
  
  velocity_key_ = plist_.get<std::string>("velocity key",
                                     Keys::getKey(domain_name,"velocity"));

  tau_e_ = plist_.get<double>("critical shear stress");
  Qe_0_ = plist_.get<double>("empirical coefficient");
  gamma_ = plist_.get<double>("specific weight of water");
  umax_ = plist_.get<double>("max current");
  xi_ = plist_.get<double>("Chezy parameter");

  double pi = boost::math::constants::pi<double>();

  lambda_ = 8./(3*pi) * (umax_/(xi_*xi_));
    
  dependencies_.insert("surface-pressure");
    
}

  
ErosionRateEvaluator ::ErosionRateEvaluator (const ErosionRateEvaluator & other) :
  EvaluatorSecondaryMonotypeCV(other),
  velocity_key_(other.velocity_key_) {

  tau_e_ = other.tau_e_;
  Qe_0_ = other.Qe_0_;
  gamma_ = other.gamma_;
  lambda_ = other.lambda_;
  umax_ = other.umax_;
  xi_ = other.xi_;
} 


Teuchos::RCP<Evaluator> ErosionRateEvaluator ::Clone() const {
  return Teuchos::rcp(new ErosionRateEvaluator (*this));
}


void ErosionRateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  const Epetra_MultiVector& vel = *S->Get<CompositeVector>(velocity_key_).ViewComponent("cell");
  Epetra_MultiVector& result_c = *result->ViewComponent("cell");
  
  for (int c=0; c<vel.MyLength(); c++){
    double tau_0 = gamma_ * lambda_ * (sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]));
    if (tau_0 > tau_e_){
      result_c[0][c] = Qe_0_*(tau_0 / tau_e_ - 1);
    }else{
      result_c[0][c] = 0.;
    }
  }

 

}

void ErosionRateEvaluator::EvaluateFieldPartialDerivative_ (const Teuchos::Ptr<State>& S,
                                                            Key wrt_key,
                                                            const Teuchos::Ptr<CompositeVector>& result) {
   AMANZI_ASSERT(0); 
}
  
  
} // namespace
