/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Suction head = \Psi( sat )

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "suction_head_evaluator.hh"

namespace Amanzi {
namespace Flow {

SuctionHeadEvaluator::SuctionHeadEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    min_val_(0.) {

  AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(sublist);
  InitializeFromPlist_();
}

SuctionHeadEvaluator::SuctionHeadEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPartition>& wrms) :
    SecondaryVariableFieldEvaluator(plist),
    wrms_(wrms),
    min_val_(0.) {
  InitializeFromPlist_();
}

SuctionHeadEvaluator::SuctionHeadEvaluator(const SuctionHeadEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrms_(other.wrms_),
    sat_key_(other.sat_key_),
    min_val_(other.min_val_) {}
 
Teuchos::RCP<FieldEvaluator>
SuctionHeadEvaluator::Clone() const {
  return Teuchos::rcp(new SuctionHeadEvaluator(*this));
}


void SuctionHeadEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("suction head key", "suction_head");
  }

  // dependencies
  Key domain_name = Keys::getDomain(my_key_);

  // -- saturation liquid
  sat_key_ = plist_.get<std::string>("saturation key",
          Keys::getKey(domain_name, "saturation_liquid"));
  dependencies_.insert(sat_key_);
  
  // cutoff above 0?
  min_val_ = plist_.get<double>("minimum suction cutoff", 0.);

}


// Special purpose EnsureCompatibility required because of surface rel perm.
void SuctionHeadEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {

    SecondaryVariableFieldEvaluator::EnsureCompatibility(S);

}


void SuctionHeadEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result->Mesh(), -1);
    wrms_->first->Verify();
  }

  // Evaluate suction.
  // -- Evaluate the model to calculate suction on cells.
  const Epetra_MultiVector& sat_c = *S->GetFieldData(sat_key_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  int ncells = res_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    int index = (*wrms_->first)[c];
    res_c[0][c] = wrms_->second[index]->suction_head(sat_c[0][c]);
  }

}


void SuctionHeadEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result->Mesh(), -1);
    wrms_->first->Verify();
  }

  if (wrt_key == sat_key_) {
    // d(psi) / dsl 
    
   
    const Epetra_MultiVector& sat_c = *S->GetFieldData(sat_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

    int ncells = res_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      int index = (*wrms_->first)[c];
      res_c[0][c] = wrms_->second[index]->d_suction_head(sat_c[0][c]);
      // AMANZI_ASSERT(res_c[0][c] >= 0.);
    }
  }
  
}



} //namespace
} //namespace
