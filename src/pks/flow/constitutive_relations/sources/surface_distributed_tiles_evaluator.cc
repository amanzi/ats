/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! Evaluates water/solute source which represent effect of distributed subsurface tiles on overland flow
/*!

Accumulated surface sources due to tile drains.

.. _surface-distributed-tiles-spec:
.. admonition:: surface-distributed-tiles-spec

   * `"number of ditches`" ``[int]`` Number of ditches, corresponding to the number of unique IDs.

   KEYS:
   - `"accumulated source`" **SUBSURFACE_DOMAIN-accumulated_source** Source to the ditch from the tile.
   - `"catchment ID`" **DOMAIN-catchments_id** ID indicating which ditch a given cell drains to.
   - `"catchment fraction`" **DOMAIN-catchments_id** 1/dL, the fraction describing the length scale used in connecting the pipe to the ditch.

*/

#include "Key.hh"
#include "surface_distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


SurfDistributedTilesRateEvaluator::SurfDistributedTilesRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  catch_id_key_ = Keys::readKey(plist, domain_, "catchment ID", "catchments_id");
  dependencies_.insert(KeyTag{ catch_id_key_, tag });

  catch_frac_key_ = Keys::readKey(plist, domain_, "catchment fraction", "catchments_frac");
  dependencies_.insert(KeyTag{ catch_frac_key_, tag });

  auto domain_subsurf = Keys::readDomainHint(plist, domain_, "surface", "domain");
  acc_sources_key_ =
    Keys::readKey(plist, domain_subsurf, "accumulated source", "accumulated_tile_sources");

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  num_ditches_ = plist.get<int>("number of ditches");
}

// Required methods from SecondaryVariableFieldEvaluator
void
SurfDistributedTilesRateEvaluator::Update_(State& S)
{
  auto key_tag = my_keys_.front();
  auto tag = key_tag.second;
  double dt = S.Get<double>("dt", tag);

  const auto& catch_id = *S.Get<CompositeVector>(catch_id_key_, tag).ViewComponent("cell", false);
  const auto& catch_frac =
    *S.Get<CompositeVector>(catch_frac_key_, tag).ViewComponent("cell", false);
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& acc_sources_vec = S.Get<Teuchos::Array<double>>(acc_sources_key_, tag);

  auto& surf_src =
    *S.GetW<CompositeVector>(key_tag.first, tag, key_tag.first).ViewComponent("cell");
  double total = 0.0;

  AmanziMesh::Entity_ID ncells = catch_id.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    if ((catch_id[0][c] > 0) && (dt > 1e-14)) {
      surf_src[0][c] = -acc_sources_vec[catch_id[0][c] - 1] * catch_frac[0][c] / (cv[0][c] * dt);
    }
  }

  total = 0.0;
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) { total += surf_src[0][c] * cv[0][c]; }
}

void
SurfDistributedTilesRateEvaluator::EnsureCompatibility(State& S)
{
  Key key = my_keys_.front().first;
  auto tag = my_keys_.front().second;

  if (!S.HasRecord(acc_sources_key_, tag)) {
    S.Require<Teuchos::Array<double>>(num_ditches_, acc_sources_key_, tag);
  }
  S.Require<CompositeVector, CompositeVectorSpace>(key, tag, key)
    .SetMesh(S.GetMesh(domain_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // Cop-out -- ensure not fully implemented for this evaluator.  FIXME --ETC
  for (const auto& dep : dependencies_) { S.RequireEvaluator(dep.first, dep.second); }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
