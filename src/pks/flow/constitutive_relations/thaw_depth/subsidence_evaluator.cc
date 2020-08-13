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
    : SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomain(my_key_);
  auto pos = domain_.find_last_of('_');
  int col_id = std::stoi(domain_.substr(pos+1, domain_.size()));
  
  std::stringstream domain_ss;
  domain_ss << "column_"<< col_id;
  bp_key_ = Keys::getKey(domain_ss.str(),"base_porosity");
  dependencies_.insert(bp_key_);

  std::stringstream domain_sf;
  domain_sf << "surface_column_"<< col_id;
  init_elev_key_ = Keys::getKey(domain_sf.str(),"initial_elevation");
  //  dependencies_.insert(init_elev_key_);
}
  

SubsidenceEvaluator::SubsidenceEvaluator(const SubsidenceEvaluator& other)
  : SecondaryVariableFieldEvaluator(other),
    bp_key_(other.bp_key_),
    init_elev_key_(other.init_elev_key_)
{}
  
Teuchos::RCP<FieldEvaluator>
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

  const auto& init_elev_c = *S->GetFieldData(init_elev_key_)
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
  bool changed = SecondaryVariableFieldEvaluator::HasFieldChanged(S,request);

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