/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
An evaluator that simply saves the initial valule of another field for later use.

Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "initial_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

InitialElevationEvaluator::InitialElevationEvaluator(Teuchos::ParameterList& plist)
  : SecondaryVariableFieldEvaluator(plist),
    updated_once_(false)
{
  Key domain_surf = Keys::getDomain(my_key_);
  Key domain_ss = Keys::readDomainHint(plist, domain_surf, "surface", "subsurface");
  bp_key_ = Keys::readKey(plist, domain_ss, "base porosity", "base_porosity");
  dependencies_.insert(bp_key_);

  // fixme -- hard-code this to checkpoint?  Or make sure this evaluates based
  // on the un-deformed initial mesh prior to re-deforming the mesh on restart?
  // I believe that it currently relies on the user to "checkpoint" this.
  // --ETC
}


Teuchos::RCP<FieldEvaluator>
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

  AMANZI_ASSERT(res_c.MyLength() == 1);
  res_c[0][0] = top_z_centroid[2];
}


// Custom EnsureCompatibility forces this to be updated once and only once.
//
// This is a bit of a hack...
bool
InitialElevationEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request)
{
  if (!updated_once_) {
    bool changed = SecondaryVariableFieldEvaluator::HasFieldChanged(S,request);
    if (changed) {
      UpdateField_(S);
      updated_once_ = true;
    }
    return changed;
  }
  return false;
}

void
InitialElevationEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  if (my_fac->Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (const auto& dep : dependencies_) {
      S->RequireField(dep)
        ->SetMesh(S->GetMesh(Keys::getDomain(dep)));
      S->RequireFieldEvaluator(dep)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
