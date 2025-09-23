/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded depth from ponded depth and subgrid parameters.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.

*/

#include "MeshAlgorithms.hh"
#include "volumetric_ponded_depth_evaluator.hh"
#include "subgrid_microtopography.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

VolumetricPondedDepthEvaluator::VolumetricPondedDepthEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // dependencies
  pd_key_ = Keys::readKey(plist_, domain, "ponded depth key", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });

  delta_max_key_ =
    Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(KeyTag{ delta_max_key_, tag });

  delta_ex_key_ = Keys::readKey(plist_, domain, "excluded volume", "excluded_volume");
  dependencies_.insert(KeyTag{ delta_ex_key_, tag });
}


void
VolumetricPondedDepthEvaluator::Evaluate_(const State& S,
                                          const std::vector<CompositeVector*>& result)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pd_v = S.GetPtr<CompositeVector>(pd_key_, tag);
  Teuchos::RCP<const CompositeVector> del_max_v = S.GetPtr<CompositeVector>(delta_max_key_, tag);
  Teuchos::RCP<const CompositeVector> del_ex_v = S.GetPtr<CompositeVector>(delta_ex_key_, tag);
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  for (const auto& comp : *result[0]) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& pd = *pd_v->ViewComponent(comp, false);
    const auto& del_max = *del_max_v->ViewComponent(internal_comp, false);
    const auto& del_ex = *del_ex_v->ViewComponent(internal_comp, false);
    auto& res = *result[0]->ViewComponent(comp, false);

    int ncomp = result[0]->size(comp, false);
    for (int i = 0; i != ncomp; ++i) {
      int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
      AMANZI_ASSERT(Microtopography::validParameters(del_max[0][ii], del_ex[0][ii]));
      res[0][i] = Microtopography::volumetricDepth(pd[0][i], del_max[0][ii], del_ex[0][ii]);
    }
  }
}


void
VolumetricPondedDepthEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pd_v = S.GetPtr<CompositeVector>(pd_key_, tag);
  Teuchos::RCP<const CompositeVector> del_max_v = S.GetPtr<CompositeVector>(delta_max_key_, tag);
  Teuchos::RCP<const CompositeVector> del_ex_v = S.GetPtr<CompositeVector>(delta_ex_key_, tag);
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  if (wrt_key == pd_key_) {
    for (const auto& comp : *result[0]) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const auto& pd = *pd_v->ViewComponent(comp, false);
      const auto& del_max = *del_max_v->ViewComponent(internal_comp, false);
      const auto& del_ex = *del_ex_v->ViewComponent(internal_comp, false);
      auto& res = *result[0]->ViewComponent(comp, false);

      int ncomp = result[0]->size(comp, false);
      for (int i = 0; i != ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        res[0][i] =
          Microtopography::dVolumetricDepth_dDepth(pd[0][i], del_max[0][ii], del_ex[0][ii]);
        res[0][i] = std::max(res[0][i], 0.001);
      }
    }

  } else {
    Errors::Message msg("VolumetricPondedDepthEvaluator: Not Implemented: no derivatives "
                        "implemented other than ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
}

void
VolumetricPondedDepthEvaluator::EnsureCompatibility_ToDeps_(State& S)
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
      if (key_tag.first == delta_max_key_ || key_tag.first == delta_ex_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*no_bf_dep_fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*dep_fac);
      }
    }
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
