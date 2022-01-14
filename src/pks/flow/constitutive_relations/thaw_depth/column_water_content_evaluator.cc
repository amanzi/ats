/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column water content evaluator gets the subsurface water content and map to surface_star 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "column_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {



ColumnWaterContentEvaluator::ColumnWaterContentEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Key dset_name = plist.get<std::string>("domain set name", "column");
  
  domain_ = Keys::getDomain(my_key_); //surface_column domain
  int col_id = Keys::getDomainSetIndex<int>(domain_);
  
  Key domain_ss = Keys::getDomainInSet(dset_name, col_id);
  
  wc_key_ = Keys::readKey(plist, domain_ss, "water content", "water_content");
  dependencies_.insert(wc_key_);
}
  

ColumnWaterContentEvaluator::ColumnWaterContentEvaluator(const ColumnWaterContentEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    wc_key_(other.wc_key_)
{}
  
Teuchos::RCP<Evaluator>
ColumnWaterContentEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnWaterContentEvaluator(*this));
}


void
ColumnWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  // search through the column and find the deepest unfrozen cell

  const auto& wc_c = *S->Get<CompositeVector>(wc_key_).ViewComponent("cell", false);

  int col_cells = wc_c.MyLength();
  double col_sum = 0;
  
  for (int i=0; i!=col_cells; ++i) {
    col_sum += wc_c[0][i];
  }
  
  res_c[0][0] = col_sum;
 
}
  
void
ColumnWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
ColumnWaterContentEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  bool changed = EvaluatorSecondaryMonotypeCV::HasFieldChanged(S,request);

  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ColumnWaterContentEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{

  AMANZI_ASSERT(my_key_ != std::string(""));
   
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT,  my_key_);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
