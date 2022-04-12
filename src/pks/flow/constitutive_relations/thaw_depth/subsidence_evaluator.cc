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
  : SecondaryVariableFieldEvaluator(plist),
    updated_once_(false)
{
  domain_ = Keys::getDomain(my_key_); //surface_column domain
  Key domain_ss = Keys::readDomainHint(plist, domain_, "surface", "subsurface");

  bp_key_ = Keys::readKey(plist, domain_ss, "base porosity", "base_porosity");
  dependencies_.insert(bp_key_);

  init_elev_key_ = Keys::readKey(plist, domain_, "initial elevation", "initial_elevation");
  dependencies_.insert(init_elev_key_);
}

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

  std::string domain_ss = Keys::getDomain(bp_key_);
  const auto& top_z_centroid = S->GetMesh(domain_ss)->face_centroid(0);
  const auto& init_elev_c = *S->GetFieldData(init_elev_key_)
    ->ViewComponent("cell", false);

  AMANZI_ASSERT(res_c.MyLength() == 1);
  AMANZI_ASSERT(init_elev_c.MyLength() == 1);
  res_c[0][0] = init_elev_c[0][0] - top_z_centroid[2];
}


void
SubsidenceEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  if (my_fac->Mesh() != Teuchos::null) {
    S->RequireField(init_elev_key_)
      ->Update(*my_fac);
    S->RequireField(bp_key_)
      ->SetMesh(S->GetMesh(Keys::getDomain(bp_key_)));

    // Recurse into the tree to propagate info to leaves.
    for (const auto& dep : dependencies_)
      S->RequireFieldEvaluator(dep)->EnsureCompatibility(S);
  }
}



} //namespace
} //namespace
