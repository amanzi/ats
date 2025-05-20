/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "Key.hh"
#include "snow_meltrate_evaluator.hh"
#include "errors.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {


SnowMeltRateEvaluator::SnowMeltRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), compatibility_checked_(false)
{
  melt_rate_ = plist.get<double>("snow melt rate [mm day^-1 C^-1]", 2.74);
  melt_rate_ *= 0.001 / 86400.; // convert mm/day to m/s

  // snow begins to melt when the air temperature is this many degrees above 0
  if (plist.isParameter("air-snow temperature difference [C]")) {
    Errors::Message msg(
      "SnowMeltRateEvaluator: parameter \"air-snow temperature difference [C]\" is no longer "
      "accepted"
      "-- instead add a new evaluator for \"snow-expected_temperature\" that is of type "
      "\"additive evaluator\" that shifts the air temperature.");
    Exceptions::amanzi_throw(msg);
  }

  Tag tag = my_keys_.front().second;
  domain_ = Keys::getDomain(my_keys_.front().first);

  exp_temp_key_ =
    Keys::readKey(plist, domain_, "expected snow temperature", "expected_temperature");
  dependencies_.insert(KeyTag{ exp_temp_key_, tag });

  snow_key_ = Keys::readKey(plist, domain_, "snow water equivalent", "water_equivalent");
  dependencies_.insert(KeyTag{ snow_key_, tag });
}

void
SnowMeltRateEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // new state!
  land_cover_ = getLandCover(S.GetModelParameters("land cover types"), { "snow_transition_depth" });
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S);
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
SnowMeltRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  auto mesh = S.GetMesh(domain_);

  const auto& exp_temp = *S.Get<CompositeVector>(exp_temp_key_, tag).ViewComponent("cell", false);
  const auto& swe = *S.Get<CompositeVector>(snow_key_, tag).ViewComponent("cell", false);
  auto& res = *result[0]->ViewComponent("cell", false);

  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh->getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : lc_ids) {
      if (exp_temp[0][c] > 273.15) {
        res[0][c] = melt_rate_ * (exp_temp[0][c] - 273.15);

        if (swe[0][c] < lc.second.snow_transition_depth) {
          res[0][c] *= std::max(0., swe[0][c] / lc.second.snow_transition_depth);
        }

      } else {
        res[0][c] = 0.0;
      }
    }
  }
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
SnowMeltRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  auto mesh = S.GetMesh(domain_);
  const auto& exp_temp = *S.Get<CompositeVector>(exp_temp_key_, tag).ViewComponent("cell", false);
  const auto& swe = *S.Get<CompositeVector>(snow_key_, tag).ViewComponent("cell", false);
  auto& res = *result[0]->ViewComponent("cell", false);

  if (wrt_key == exp_temp_key_) {
    for (const auto& lc : land_cover_) {
      auto lc_ids = mesh->getSetEntities(
        lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      for (auto c : lc_ids) {
        if (exp_temp[0][c] > 273.15) {
          res[0][c] = melt_rate_;
          if (swe[0][c] < lc.second.snow_transition_depth) {
            res[0][c] *= std::max(0., swe[0][c] / lc.second.snow_transition_depth);
          }
        } else {
          res[0][c] = 0.0;
        }
      }
    }

  } else if (wrt_key == snow_key_) {
    for (const auto& lc : land_cover_) {
      auto lc_ids = mesh->getSetEntities(
        lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      for (auto c : lc_ids) {
        if (swe[0][c] < lc.second.snow_transition_depth && exp_temp[0][c] > 273.15) {
          res[0][c] = melt_rate_ * (exp_temp[0][c] - 273.15) / lc.second.snow_transition_depth;
        } else {
          res[0][c] = 0.0;
        }
      }
    }
  }
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
