/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "organic_matter_evaluator.hh"
#include "boost/math/constants/constants.hpp"

namespace Amanzi {

OrganicMatterRateEvaluator :: OrganicMatterRateEvaluator(Teuchos::ParameterList& plist) :
  EvaluatorSecondaryMonotypeCV(plist) {

  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  
  biomass_key_ = plist_.get<std::string>("biomass key", Keys::getKey(domain_name,"biomass"));

  Bmax_ = plist_.get<double>("maximum biomass");
  Q_db0_ = plist_.get<double>("empirical coefficient");
 
   
  dependencies_.insert(KeyTag{biomass_key_, tag});
    
}

  
OrganicMatterRateEvaluator ::OrganicMatterRateEvaluator (const OrganicMatterRateEvaluator & other) :
  EvaluatorSecondaryMonotypeCV(other) {

  biomass_key_ = other.biomass_key_;
  Bmax_ = other.Bmax_;
  Q_db0_ = other.Q_db0_;

} 


Teuchos::RCP<Evaluator> OrganicMatterRateEvaluator ::Clone() const {
  return Teuchos::rcp(new OrganicMatterRateEvaluator (*this));
}


void OrganicMatterRateEvaluator::Evaluate_(const State& S,
                                           const std::vector<CompositeVector*>& result){

  Tag tag = my_keys_.front().second;
  
  const Epetra_MultiVector& bio = *S.GetPtr<CompositeVector>(biomass_key_, tag)->ViewComponent("cell");
  Epetra_MultiVector& result_c = *result[0]->ViewComponent("cell");

  result_c.PutScalar(0.);
  
  for (int c=0; c<result_c.MyLength(); c++){
    for (int j=0; j<bio.NumVectors(); j++){
      result_c[0][c] +=  Q_db0_*bio[j][c]/Bmax_;
    }
  }

 

}

void OrganicMatterRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                            const Key& wrt_key, const Tag& wrt_tag,
                                                            const std::vector<CompositeVector*>& result){
                                 
   AMANZI_ASSERT(0); 
}
  
  
} // namespace
