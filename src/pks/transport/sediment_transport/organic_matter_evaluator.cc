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
  SecondaryVariableFieldEvaluator(plist) {

  Key domain_name = "surface";
  
  biomass_key_ = plist_.get<std::string>("biomass key", Keys::getKey(domain_name,"biomass"));

  Bmax_ = plist_.get<double>("maximum biomass");
  Q_db0_ = plist_.get<double>("empirical coefficient");
 
   
  dependencies_.insert(biomass_key_);
    
}

  
OrganicMatterRateEvaluator ::OrganicMatterRateEvaluator (const OrganicMatterRateEvaluator & other) :
  SecondaryVariableFieldEvaluator(other) {

  biomass_key_ = other.biomass_key_;
  Bmax_ = other.Bmax_;
  Q_db0_ = other.Q_db0_;

} 


Teuchos::RCP<FieldEvaluator> OrganicMatterRateEvaluator ::Clone() const {
  return Teuchos::rcp(new OrganicMatterRateEvaluator (*this));
}


void OrganicMatterRateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  const Epetra_MultiVector& bio = *S->Get<CompositeVector>(biomass_key_).ViewComponent("cell");
  Epetra_MultiVector& result_c = *result->ViewComponent("cell");

  result_c.PutScalar(0.);
  
  for (int c=0; c<result_c.MyLength(); c++){
    for (int j=0; j<bio.NumVectors(); j++){
      result_c[0][c] +=  Q_db0_*bio[j][c]/Bmax_;
    }
  }

 

}

void OrganicMatterRateEvaluator::EvaluateFieldPartialDerivative_ (const Teuchos::Ptr<State>& S,
                                                            Key wrt_key,
                                                            const Teuchos::Ptr<CompositeVector>& result) {
   AMANZI_ASSERT(0); 
}
  
  
} // namespace
