/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: 
*/

/*!

Requires the following dependencies:


*/

#include "Key.hh"
#include "surface_distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


SurfDistributedTilesRateEvaluator::SurfDistributedTilesRateEvaluator(Teuchos::ParameterList& plist) :
  EvaluatorSecondaryMonotypeCV(plist),
  compatibility_checked_(false)
{

  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  surface_marks_key_ = Keys::readKey(plist, domain_, "catchments_id", "catchments_id");
  surf_len_key_ = Keys::readKey(plist, domain_, "catchments_frac", "catchments_frac");
  dist_sources_key_ = plist.get<std::string>("accumulated source key", "subdomain_sources");
  Key update_key = plist.get<std::string>("update key", "surface-pressure");
  
  num_ditches_ = plist.get<int>("number of ditches");
  implicit_ = plist.get<bool>("implicit drainage", true);
  
  dependencies_.insert(KeyTag{surface_marks_key_, tag});
  dependencies_.insert(KeyTag{surf_len_key_, tag});
  dependencies_.insert(KeyTag{update_key, tag});
 
  // times_.resize(2);
  // times_[0] = -1.0; times_[1] = -1.0;
  
}

// Required methods from SecondaryVariableFieldEvaluator
void
SurfDistributedTilesRateEvaluator::Evaluate_(const State& S,
                                                  const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second; 
  double dt = S.Get<double>("dt", tag); 
    
  // Teuchos::RCP<Field> src_field =  S->GetField(dist_sources_key_, "state");
  // if (!implicit_){
  //   if (abs(t0 - times_[0]) > 1e-12) {
  //     src_field->SwitchCopies(Key("default"), Key("prev_timestep"));
  //     times_[0] = t0;
  //   }
  // }
      
  const auto& surf_marks = *S.GetPtr<CompositeVector>(surface_marks_key_, tag)->ViewComponent("cell", false);
  const auto& len_frac = *S.GetPtr<CompositeVector>(surf_len_key_, tag)->ViewComponent("cell", false);
  const auto& cv =
    *S.GetPtr<CompositeVector>(Keys::getKey(domain_,"cell_volume"), tag)->ViewComponent("cell",false);  
  auto& surf_src = *result[0]->ViewComponent("cell",false);
  
  auto dist_src_vec = S.GetPtr<Epetra_Vector>(dist_sources_key_, tag);
                                            
  // if (implicit_){
  //   src_vec = src_field->GetConstantVectorData();
  // }else{
  //   src_vec = src_field->GetCopy("prev_timestep")->GetConstantVectorData();
  // }
  
  //std::cout << "Surface source\n"<< *src_vec<<"\n";

  double total = 0.0;
  double frac_len = 0.0;
  
  AmanziMesh::Entity_ID ncells = surf_marks.MyLength();
  for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
    if ((surf_marks[0][c] > 0) && (dt>1e-14)) {
      surf_src[0][c] = -(*dist_src_vec)[surf_marks[0][c] - 1]*len_frac[0][c]/(cv[0][c]*dt) ;
      //surf_src[0][c] = 0.0;
      if (surf_marks[0][c]==1) frac_len += len_frac[0][c];
        
      total += -(*dist_src_vec)[surf_marks[0][c]-1]*len_frac[0][c];
    }
  }
  
  // std::cout<< "Total source: "<<total<<" "<<frac_len<<"\n";
  // total = 0.0;
  // for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
  //   total += surf_src[0][c]*cv[0][c];
  // }
  // std::cout<<"Total source: field "<<my_key_<<" "<<total<<"\n";
}

// Required methods from SecondaryVariableFieldEvaluator
void
SurfDistributedTilesRateEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{

  result[0]->PutScalar(0.0);
  
}

void
SurfDistributedTilesRateEvaluator::EnsureCompatibility_Structure_(State& S) {

  auto tag = my_keys_.front().second;
  if (!S.HasRecord(dist_sources_key_, tag)){    
  //   S->RequireConstantVector(dist_sources_key_, num_ditches_);
  //   Teuchos::RCP<Field> field =  S->GetField(dist_sources_key_, "state");
  //   Teuchos::RCP<Field_ConstantVector> cvfield =
  //     Teuchos::rcp_dynamic_cast<Field_ConstantVector>(field, true);
  //   cvfield->CreateData();

  //   field->GetConstantVectorData()->PutScalar(0.0);
  //   field->set_initialized();
    
  //   Key prev_step("prev_timestep");
  //   field->RequireCopy(prev_step);
  }


  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Structure_(S);

  
}

void
SurfDistributedTilesRateEvaluator::InitializeFromPlist_() {

}



} //namespace
} //namespace
} //namespace

