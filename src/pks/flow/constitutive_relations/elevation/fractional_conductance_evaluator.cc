/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "MeshAlgorithms.hh"
#include "fractional_conductance_evaluator.hh"
#include "subgrid_microtopography.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace FlowRelations {

FractionalConductanceEvaluator::FractionalConductanceEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  vpd_key_ = Keys::readKey(plist_, domain, "volumetric ponded depth", "volumetric_ponded_depth");
  dependencies_.insert(KeyTag{ vpd_key_, tag });

  mobile_depth_key_ = Keys::readKey(plist_, domain, "mobile depth", "mobile_depth");
  dependencies_.insert(KeyTag{ mobile_depth_key_, tag });

  delta_max_key_ =
    Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(KeyTag{ delta_max_key_, tag });

  delta_ex_key_ = Keys::readKey(plist_, domain, "excluded volume", "excluded_volume");
  dependencies_.insert(KeyTag{ delta_ex_key_, tag });

  depr_depth_key_ = Keys::readKey(plist_, domain, "depression depth", "depression_depth");
  dependencies_.insert(KeyTag{ depr_depth_key_, tag });
}


Teuchos::RCP<Evaluator>
FractionalConductanceEvaluator::Clone() const
{
  return Teuchos::rcp(new FractionalConductanceEvaluator(*this));
}


void
FractionalConductanceEvaluator::Evaluate_(const State& S,
                                          const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> vpd_v = S.GetPtr<CompositeVector>(vpd_key_, tag);
  Teuchos::RCP<const CompositeVector> del_max_v = S.GetPtr<CompositeVector>(delta_max_key_, tag);
  Teuchos::RCP<const CompositeVector> del_ex_v = S.GetPtr<CompositeVector>(delta_ex_key_, tag);
  Teuchos::RCP<const CompositeVector> depr_depth_v =
    S.GetPtr<CompositeVector>(depr_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> mobile_depth_v =
    S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  const auto& mesh = *result[0]->Mesh();

  for (const auto& comp : *result[0]) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& vpd = *vpd_v->ViewComponent(comp, false);
    const auto& mobile_depth = *mobile_depth_v->ViewComponent(comp, false);

    const auto& del_max = *del_max_v->ViewComponent(internal_comp, false);
    const auto& del_ex = *del_ex_v->ViewComponent(internal_comp, false);
    const auto& depr_depth = *depr_depth_v->ViewComponent(internal_comp, false);
    auto& res = *result[0]->ViewComponent(comp, false);

    int ncomp = result[0]->size(comp, false);
    for (int i = 0; i != ncomp; ++i) {
      int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
      if (mobile_depth[0][i] <= 0.0) {
        res[0][i] = 0;
      } else {
        double vol_depr_depth =
          Microtopography::volumetricDepth(depr_depth[0][ii], del_max[0][ii], del_ex[0][ii]);
        res[0][i] = (vpd[0][i] - vol_depr_depth) / (mobile_depth[0][i]);
      }
    }
  }
}


void
FractionalConductanceEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> vpd_v = S.GetPtr<CompositeVector>(vpd_key_, tag);
  Teuchos::RCP<const CompositeVector> del_max_v = S.GetPtr<CompositeVector>(delta_max_key_, tag);
  Teuchos::RCP<const CompositeVector> del_ex_v = S.GetPtr<CompositeVector>(delta_ex_key_, tag);
  Teuchos::RCP<const CompositeVector> depr_depth_v =
    S.GetPtr<CompositeVector>(depr_depth_key_, tag);
  Teuchos::RCP<const CompositeVector> mobile_depth_v =
    S.GetPtr<CompositeVector>(mobile_depth_key_, tag);
  const auto& mesh = *result[0]->Mesh();

  if (wrt_key == mobile_depth_key_) {
    for (const auto& comp : *result[0]) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const auto& vpd = *vpd_v->ViewComponent(comp, false);
      const auto& mobile_depth = *mobile_depth_v->ViewComponent(comp, false);

      const auto& del_max = *del_max_v->ViewComponent(internal_comp, false);
      const auto& del_ex = *del_ex_v->ViewComponent(internal_comp, false);
      const auto& depr_depth = *depr_depth_v->ViewComponent(internal_comp, false);
      auto& res = *result[0]->ViewComponent(comp, false);

      int ncomp = result[0]->size(comp, false);
      for (int i = 0; i != ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        if (mobile_depth[0][i] <= 0.0) {
          res[0][i] = 0;
        } else {
          double vol_depr_depth =
            Microtopography::volumetricDepth(depr_depth[0][ii], del_max[0][ii], del_ex[0][ii]);
          res[0][i] = -(vpd[0][i] - vol_depr_depth) * std::pow(mobile_depth[0][ii], -2.);
        }
      }
    }
  } else if (wrt_key == vpd_key_) {
    for (const auto& comp : *result[0]) {
      const auto& mobile_depth = *mobile_depth_v->ViewComponent(comp, false);
      auto& res = *result[0]->ViewComponent(comp, false);

      int ncomp = result[0]->size(comp, false);
      for (int i = 0; i != ncomp; ++i) {
        if (mobile_depth[0][i] <= 0.0) {
          res[0][i] = 0;
        } else {
          res[0][i] = 1. / mobile_depth[0][i];
        }
      }
    }
  } else {
    Errors::Message msg("FractionalConductanceEvaluator: Not Implemented: no derivatives "
                        "implemented other than mobile depth and volumentric ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
}


void
FractionalConductanceEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag my_tag = my_keys_.front().second;
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key, my_tag);

  if (my_fac.Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    auto dep_fac = Teuchos::rcp(new CompositeVectorSpace(my_fac));
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
      if (key_tag.first == depr_depth_key_ || key_tag.first == delta_ex_key_ ||
          key_tag.first == delta_max_key_) {
        auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second);
        fac.Update(*no_bf_dep_fac);
      } else {
        auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second);
        fac.Update(*dep_fac);
      }
    }
  }
}


} // namespace FlowRelations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
