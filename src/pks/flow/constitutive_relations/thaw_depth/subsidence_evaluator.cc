/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the subsurface temperature and computes the thaw depth 
  over time.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "subsidence_evaluator.hh"

namespace Amanzi {
namespace Flow {



SubsidenceEvaluator::SubsidenceEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  Key dset_name = plist.get<std::string>("domain set name", "column");
  Key surf_dset_name = plist.get<std::string>("surface domain set name", "surface_column");
  
  domain_ = Keys::getDomain(my_key_); //surface_column domain
  int col_id = Keys::getDomainSetIndex<int>(domain_);
  
  Key domain_ss = Keys::getDomainInSet(dset_name, col_id);
  Key domain_sf = Keys::getDomainInSet(surf_dset_name, col_id);

  bp_key_ = Keys::readKey(plist, domain_ss, "base porosity", "base_porosity");
  dependencies_.insert(bp_key_);

  init_elev_key_ = Keys::readKey(plist, domain_sf, "initial elevation", "initial_elevation");
}
  

SubsidenceEvaluator::SubsidenceEvaluator(const SubsidenceEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    bp_key_(other.bp_key_),
    init_elev_key_(other.init_elev_key_)
{}
  
Teuchos::RCP<Evaluator>
SubsidenceEvaluator::Clone() const
{
  return Teuchos::rcp(new SubsidenceEvaluator(*this));
}


void
SubsidenceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{ 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  
  // search through the column and find the deepest unfrozen cell

  std::string domain_ss = Keys::getDomain(bp_key_);
  const auto& top_z_centroid = S->GetMesh(domain_ss)->face_centroid(0);

  const auto& init_elev_c = *S->GetPtr<CompositeVector>(init_elev_key_)
    ->ViewComponent("cell", false);
      
  res_c[0][0] = init_elev_c[0][0] - top_z_centroid[2];

}
  
void
SubsidenceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{}

 
// Custom EnsureCompatibility forces this to be updated once.
bool
SubsidenceEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
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
SubsidenceEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
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
