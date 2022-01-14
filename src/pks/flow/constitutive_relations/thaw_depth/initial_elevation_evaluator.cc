/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes the thaw depth 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "initial_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {



InitialElevationEvaluator::InitialElevationEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Key dset_name = plist.get<std::string>("domain set name", "column");
  
  domain_ = Keys::getDomain(my_key_); //surface_column domain
  int col_id = Keys::getDomainSetIndex<int>(domain_);
  
  Key domain_ss = Keys::getDomainInSet(dset_name, col_id);

  bp_key_ = Keys::readKey(plist, domain_ss, "base porosity", "base_porosity");

}
  

InitialElevationEvaluator::InitialElevationEvaluator(const InitialElevationEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    bp_key_(other.bp_key_)
{}
  
Teuchos::RCP<Evaluator>
InitialElevationEvaluator::Clone() const
{
  return Teuchos::rcp(new InitialElevationEvaluator(*this));
}


void
InitialElevationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  // search through the column and find the deepest unfrozen cell

  std::string domain_ss = Keys::getDomain(bp_key_);
  const auto& top_z_centroid = S->GetMesh(domain_ss)->face_centroid(0);
  
  res_c[0][0] = std::max(res_c[0][0],top_z_centroid[2]);
  
}
  
void
InitialElevationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
InitialElevationEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
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
InitialElevationEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
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
