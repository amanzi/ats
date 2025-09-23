/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "MeshAlgorithms.hh"
#include "subgrid_mobile_depth_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace FlowRelations {

SubgridMobileDepthEvaluator::SubgridMobileDepthEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  depth_key_ = Keys::readKey(plist_, domain_name, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ depth_key_, tag });

  depr_depth_key_ = Keys::readKey(plist_, domain_name, "depression depth key", "depression_depth");
  dependencies_.insert(KeyTag{ depr_depth_key_, tag });
}


Teuchos::RCP<Evaluator>
SubgridMobileDepthEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridMobileDepthEvaluator(*this));
}


void
SubgridMobileDepthEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  auto depr_depth_v = S.GetPtr<CompositeVector>(depr_depth_key_, tag);
  auto depth_v = S.GetPtr<CompositeVector>(depth_key_, tag);
  const auto& mesh = *result[0]->Mesh();

  for (const auto& comp : *result[0]) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& depth = *depth_v->ViewComponent(comp, false);
    const auto& depr_depth = *depr_depth_v->ViewComponent(internal_comp, false);
    auto& res = *result[0]->ViewComponent(comp, false);

    int ncomp = result[0]->size(comp, false);
    for (int i = 0; i != ncomp; ++i) {
      int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
      res[0][i] = std::max(0., depth[0][i] - depr_depth[0][ii]);
    }
  }
}


void
SubgridMobileDepthEvaluator::EvaluatePartialDerivative_(const State& S,
                                                        const Key& wrt_key,
                                                        const Tag& wrt_tag,
                                                        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  auto depr_depth_v = S.GetPtr<CompositeVector>(depr_depth_key_, tag);
  auto depth_v = S.GetPtr<CompositeVector>(depth_key_, tag);
  const auto& mesh = *result[0]->Mesh();

  if (wrt_key == depth_key_) {
    result[0]->PutScalar(1.);
  } else {
    result[0]->PutScalar(-1.);
  }
}


void
SubgridMobileDepthEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // Ensure my field exists.  Requirements should be already set.
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                        my_keys_.front().second);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac = Teuchos::rcp(new CompositeVectorSpace(my_fac));
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
    for (const auto& key_tag : dependencies_) {
      if (key_tag.first == depr_depth_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*no_bf_dep_fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*dep_fac);
      }
    }
  }
}


} // namespace FlowRelations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
