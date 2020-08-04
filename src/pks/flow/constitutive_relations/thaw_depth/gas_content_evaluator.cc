/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The gas evaluator gets the subsurface temperature and gas saturation and computes gas content 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "gas_content_evaluator.hh"

namespace Amanzi {
namespace Flow {



GasContentEvaluator::GasContentEvaluator(Teuchos::ParameterList& plist)
    : SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomain(my_key_);
  auto pos = domain_.find_last_of('_');
  int col_id = std::stoi(domain_.substr(pos+1, domain_.size()));
  
  std::stringstream domain_ss;
  domain_ss << "column_"<< col_id;
  temp_key_ = Keys::getKey(domain_ss.str(),"temperature");
  dependencies_.insert(temp_key_);
  
  cv_key_ = Keys::getKey(domain_ss.str(),"cell_volume");
  dependencies_.insert(cv_key_);
  
  sat_key_ = Keys::getKey(domain_ss.str(),"saturation_gas");
  dependencies_.insert(sat_key_);
  
  trans_width_ =  plist_.get<double>("transition width [K]", 0.2);
}
  

GasContentEvaluator::GasContentEvaluator(const GasContentEvaluator& other)
  : SecondaryVariableFieldEvaluator(other),
    temp_key_(other.temp_key_),
    cv_key_(other.cv_key_),
    sat_key_(other.sat_key_)
{}
  
Teuchos::RCP<FieldEvaluator>
GasContentEvaluator::Clone() const
{
  return Teuchos::rcp(new GasContentEvaluator(*this));
}


void
GasContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  double trans_temp = 273.15 + 0.5*trans_width_;

  // search through the column and find the deepest unfrozen cell

  std::string domain_ss = Keys::getDomain(temp_key_);
  const auto& top_z_centroid = S->GetMesh(domain_ss)->face_centroid(0);
  AmanziGeometry::Point z_centroid(top_z_centroid);

  const auto& temp_c = *S->GetFieldData(temp_key_)->ViewComponent("cell", false);
  const auto& cv_c = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);
  const auto& sat_c = *S->GetFieldData(sat_key_)->ViewComponent("cell", false);

  int col_cells = temp_c.MyLength();
  double col_sum = 0;
  double cv_sum = 0;
  for (int i=0; i!=col_cells; ++i) {
    // this hard codes in the transition width to 0.2 K
    if (sat_c[0][i]>0) {
      col_sum += cv_c[0][i] * sat_c[0][i];
      cv_sum += cv_c[0][i];
    }
  }
  
  res_c[0][0] = cv_sum > 0.0 ? col_sum/cv_sum : 0;
 
}
  
void
GasContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
GasContentEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  bool changed = SecondaryVariableFieldEvaluator::HasFieldChanged(S,request);

  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
GasContentEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{

  AMANZI_ASSERT(my_key_ != std::string(""));
   
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
