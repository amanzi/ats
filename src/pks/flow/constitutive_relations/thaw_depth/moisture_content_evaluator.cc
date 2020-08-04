/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes moisture content 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "moisture_content_evaluator.hh"

namespace Amanzi {
namespace Flow {



MoistureContentEvaluator::MoistureContentEvaluator(Teuchos::ParameterList& plist)
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
  
  sat_key_ = Keys::getKey(domain_ss.str(),"saturation_liquid");
  dependencies_.insert(sat_key_);
  
  trans_width_ =  plist_.get<double>("transition width [K]", 0.2);

  volumetric_wc_ =  plist_.get<bool>("volumetric water content", false);

  if (volumetric_wc_) {
    por_key_ = Keys::getKey(domain_ss.str(),"base_porosity");
    dependencies_.insert(por_key_);
  }
}
  

MoistureContentEvaluator::MoistureContentEvaluator(const MoistureContentEvaluator& other)
  : SecondaryVariableFieldEvaluator(other),
    temp_key_(other.temp_key_),
    cv_key_(other.cv_key_),
    sat_key_(other.sat_key_),
    por_key_(other.por_key_),
    trans_width_(other.trans_width_),
    volumetric_wc_(other.volumetric_wc_)
    
{}
  
Teuchos::RCP<FieldEvaluator>
MoistureContentEvaluator::Clone() const
{
  return Teuchos::rcp(new MoistureContentEvaluator(*this));
}


void
MoistureContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  double trans_temp = 273.15 + 0.5*trans_width_;

  // search through the column and find the deepest unfrozen cell

  std::string domain_ss = Keys::getDomain(temp_key_);

  const auto& temp_c = *S->GetFieldData(temp_key_)->ViewComponent("cell", false);
  const auto& cv_c = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);
  const auto& sat_c = *S->GetFieldData(sat_key_)->ViewComponent("cell", false);

  int col_cells = temp_c.MyLength();
  double col_sum = 0;
  double cv_sum = 0;
  if (!volumetric_wc_) {
    for (int i=0; i!=col_cells; ++i) {
      if (temp_c[0][i] >= trans_temp) {
        col_sum += cv_c[0][i] * sat_c[0][i];
        cv_sum += cv_c[0][i];
      }
    }
  }
  else {
    const auto& por_c = *S->GetFieldData(por_key_)->ViewComponent("cell", false);
    for (int i=0; i!=col_cells; ++i) {
      if (temp_c[0][i] >= trans_temp) {
        col_sum += por_c[0][i] * sat_c[0][i];
        cv_sum += 1.0;
      }
    }
  }
  
  res_c[0][0] = cv_sum > 0.0 ? col_sum/cv_sum : 0;
 
}
  
void
MoistureContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
MoistureContentEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
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
MoistureContentEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
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
