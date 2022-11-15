/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*

  Evaluates water/solute source which represent effect of distributed subsurface tiles on overland flow

  License: see $ATS_DIR/COPYRIGHT
  Author:
*/

/*!

Requires the following dependencies:


*/

#include "Key.hh"
#include "surface_distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


SurfDistributedTilesRateEvaluator::SurfDistributedTilesRateEvaluator(Teuchos::ParameterList& plist) :
  EvaluatorSecondaryMonotypeCV(plist),
  compatibility_checked_(false)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  surface_marks_key_ = Keys::readKey(plist, domain_, "catchments_id", "catchments_id");
  surf_len_key_ = Keys::readKey(plist, domain_, "catchments_frac", "catchments_frac");
  dist_sources_key_ = plist.get<std::string>("accumulated source key", "subdomain_sources");
  Key update_key = plist.get<std::string>("update key", "water_source");

  num_ditches_ = plist.get<int>("number of ditches");
  implicit_ = plist.get<bool>("implicit drainage", true);

  dependencies_.insert(KeyTag{surface_marks_key_, tag});
  dependencies_.insert(KeyTag{surf_len_key_, tag});
  dependencies_.insert(KeyTag{update_key, tag});
}

// Required methods from SecondaryVariableFieldEvaluator
void
SurfDistributedTilesRateEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  double dt = S.Get<double>("dt", tag);

  const auto& surf_marks = *S.Get<CompositeVector>(surface_marks_key_, tag).ViewComponent("cell", false);
  const auto& len_frac = *S.Get<CompositeVector>(surf_len_key_, tag).ViewComponent("cell", false);
  const auto& cv = *S.Get<CompositeVector>(Keys::getKey(domain_,"cell_volume"), tag).ViewComponent("cell",false);
  const auto& dist_src_vec = S.Get<Teuchos::Array<double>>(dist_sources_key_, tag);
  auto& surf_src = *result[0]->ViewComponent("cell",false);

  double total = 0.0;

  AmanziMesh::Entity_ID ncells = surf_marks.MyLength();
  for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
    if ((surf_marks[0][c] > 0) && (dt>1e-14)) {
      surf_src[0][c] = -dist_src_vec[surf_marks[0][c] - 1] * len_frac[0][c] / (cv[0][c] * dt) ;
    }
  }

  total = 0.0;
  for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
    total += surf_src[0][c] * cv[0][c];
  }
}

// Required methods from SecondaryVariableFieldEvaluator
void
SurfDistributedTilesRateEvaluator::EvaluatePartialDerivative_(const State& S,
    const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(0.0);
}

void
SurfDistributedTilesRateEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  auto tag = my_keys_.front().second;
  if (!S.HasRecord(dist_sources_key_, tag)) {
    S.Require<Teuchos::Array<double>>(num_ditches_, dist_sources_key_, tag);
  }
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S);
}

} //namespace
} //namespace
} //namespace

