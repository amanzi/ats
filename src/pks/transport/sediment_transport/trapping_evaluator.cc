/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors:
*/

#include "trapping_evaluator.hh"

namespace Amanzi {

TrappingRateEvaluator :: TrappingRateEvaluator(Teuchos::ParameterList& plist) :
  EvaluatorSecondaryMonotypeCV(plist) {

  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  
  velocity_key_ = plist_.get<std::string>("velocity key",
                                          Keys::getKey(domain_name,"velocity"));
  
  sediment_key_ =  plist_.get<std::string>("sediment key",
                                           Keys::getKey(domain_name,"sediment"));

  ponded_depth_key_ =  plist_.get<std::string>("ponded depth key",
                                           Keys::getKey(domain_name,"ponded_depth"));

  biomass_key_ = Keys::getKey(domain_name, "biomass");
  
  visc_ = plist_.get<double>("kinematic viscosity");
  d_p_ = plist_.get<double>("particle diameter");
  alpha_ = plist_.get<double>("alpha");
  beta_ = plist_.get<double>("beta");
  gamma_ = plist_.get<double>("gamma");
  sediment_density_ = plist_.get<double>("sediment density [kg m^-3]");  

  dependencies_.insert(KeyTag{"surface-pressure", tag});
  dependencies_.insert(KeyTag{sediment_key_, tag});
  dependencies_.insert(KeyTag{ponded_depth_key_, tag});
  dependencies_.insert(KeyTag{biomass_key_, tag});
    
}

  
TrappingRateEvaluator ::TrappingRateEvaluator (const TrappingRateEvaluator & other) :
  EvaluatorSecondaryMonotypeCV(other),
  velocity_key_(other.velocity_key_),
  sediment_key_(other.sediment_key_),
  ponded_depth_key_(other.ponded_depth_key_)
{
  visc_ = other.visc_;
  d_p_ = other.d_p_;
  alpha_ = other.alpha_;
  gamma_ = other.gamma_;
  beta_  = other.beta_;
  sediment_density_ = other.sediment_density_;
} 


Teuchos::RCP<Evaluator> TrappingRateEvaluator ::Clone() const {
  return Teuchos::rcp(new TrappingRateEvaluator (*this));
}


void TrappingRateEvaluator::Evaluate_(const State& S,
                                      const std::vector<CompositeVector*>& result){

  Tag tag = my_keys_.front().second;
  
  const Epetra_MultiVector& vel = *S.GetPtr<CompositeVector>(velocity_key_, tag)->ViewComponent("cell");
  const Epetra_MultiVector& tcc = *S.GetPtr<CompositeVector>(sediment_key_, tag)->ViewComponent("cell");
  const Epetra_MultiVector& depth = *S.GetPtr<CompositeVector>(ponded_depth_key_, tag)->ViewComponent("cell");
  const Epetra_MultiVector& bio_n = *S.GetPtr<CompositeVector>("surface-stem_density", tag)->ViewComponent("cell");
  const Epetra_MultiVector& bio_d = *S.GetPtr<CompositeVector>("surface-stem_diameter", tag)->ViewComponent("cell");
  const Epetra_MultiVector& bio_h = *S.GetPtr<CompositeVector>("surface-stem_height", tag)->ViewComponent("cell");
  Epetra_MultiVector& result_c = *result[0]->ViewComponent("cell");

  result_c.PutScalar(0.);
  
  for (int c=0; c<result_c.MyLength(); c++){

    for (int j=0; j<bio_n.NumVectors(); j++){

      double h_s = bio_h[j][c];
      double d_s = bio_d[j][c];
      double n_s = bio_n[j][c];
    
      double u_abs = sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]);
      double eps;
      if (d_s > 1e-12){
        eps = alpha_ * std::pow(u_abs*d_s/visc_, beta_) * std::pow(d_p_/d_s, gamma_);
      }else{
        eps = 0.;
      }      

      result_c[0][c] += sediment_density_ * tcc[0][c] * u_abs * eps * d_s * n_s * std::min(depth[0][c], h_s);
    }
  }
      

}

void TrappingRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key, const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& result){
   AMANZI_ASSERT(0); 
}
  
  
} // namespace
