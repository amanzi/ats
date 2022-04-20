/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Mesh_Algorithms.hh"
#include "subgrid_mobile_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

SubgridMobileDepthEvaluator::SubgridMobileDepthEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);

  depth_key_ = Keys::readKey(plist_, domain, "ponded depth", "ponded_depth");
  dependencies_.insert(depth_key_);

  depr_depth_key_ = plist_.get<std::string>("depression depth key", Keys::getKey(domain,"depression_depth"));
  dependencies_.insert(depr_depth_key_);
}


Teuchos::RCP<FieldEvaluator>
SubgridMobileDepthEvaluator::Clone() const {
  return Teuchos::rcp(new SubgridMobileDepthEvaluator(*this));
}


void SubgridMobileDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> depr_depth_v = S->GetFieldData(depr_depth_key_);
  Teuchos::RCP<const CompositeVector> depth_v = S->GetFieldData(depth_key_);
  const auto& mesh = *result->Mesh();

  for (const auto& comp : *result) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& depth = *depth_v->ViewComponent(comp,false);
    const auto& depr_depth = *depr_depth_v->ViewComponent(internal_comp,false);
    auto& res = *result->ViewComponent(comp,false);

    int ncomp = result->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
      res[0][i] = std::max(0., depth[0][i] - depr_depth[0][ii]);
    }
  }
}


void
SubgridMobileDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> depr_depth_v = S->GetFieldData(depr_depth_key_);
  Teuchos::RCP<const CompositeVector> depth_v = S->GetFieldData(depth_key_);
  const auto& mesh = *result->Mesh();

  if (wrt_key == depth_key_) {
    result->PutScalar(1.);
  } else {
    result->PutScalar(-1.);
  }
}


void SubgridMobileDepthEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_key_ != std::string(""));
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac->Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac =
        Teuchos::rcp(new CompositeVectorSpace(*my_fac));
    dep_fac->SetOwned(false);

    Teuchos::RCP<CompositeVectorSpace> no_bf_dep_fac;
    if (dep_fac->HasComponent("boundary_face")) {
      no_bf_dep_fac = Teuchos::rcp(new CompositeVectorSpace());
      no_bf_dep_fac->SetMesh(dep_fac->Mesh())
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    } else {
      no_bf_dep_fac = dep_fac;
    }

    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& key : dependencies_) {
      if (key == my_key_) {
        Errors::Message msg;
        msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
        Exceptions::amanzi_throw(msg);
      }

      if (key == depr_depth_key_) {
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
        fac->Update(*no_bf_dep_fac);
      } else {
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
        fac->Update(*dep_fac);
      }
    }

    // Recurse into the tree to propagate info to leaves.
    for (const auto& key : dependencies_) {
      S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
    }
  }
}


} //namespace
} //namespace
} //namespace
