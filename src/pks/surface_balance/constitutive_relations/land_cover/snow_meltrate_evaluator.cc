/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "Key.hh"
#include "snow_meltrate_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {


SnowMeltRateEvaluator::SnowMeltRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), compatibility_checked_(false)
{
  melt_rate_ = plist.get<double>("snow melt rate [mm day^-1 C^-1]", 2.74) * 0.001 /
               86400.; // convert mm/day to m/s
  snow_temp_shift_ =
    plist.get<double>("air-snow temperature difference [C]",
                      2.0); // snow is typically a few degrees colder than air at melt time

  Tag tag = my_keys_.front().second;
  domain_ = Keys::getDomain(my_keys_.front().first);
  domain_surf_ = Keys::readDomainHint(plist_, domain_, "snow", "surface");

  temp_key_ = Keys::readKey(plist, domain_surf_, "air temperature", "air_temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  snow_key_ = Keys::readKey(plist, domain_, "snow water equivalent", "water_equivalent");
  dependencies_.insert(KeyTag{ snow_key_, tag });
}

void
SnowMeltRateEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // new state!
  land_cover_ = getLandCover(S.ICList().sublist("land cover types"), { "snow_transition_depth" });
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S);
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
SnowMeltRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  auto mesh = S.GetMesh(domain_);

  const auto& air_temp = *S.Get<CompositeVector>(temp_key_, tag).ViewComponent("cell", false);
  const auto& swe = *S.Get<CompositeVector>(snow_key_, tag).ViewComponent("cell", false);
  auto& res = *result[0]->ViewComponent("cell", false);

  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh->getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : lc_ids) {
      if (air_temp[0][c] - snow_temp_shift_ > 273.15) {
        res[0][c] = melt_rate_ * (air_temp[0][c] - snow_temp_shift_ - 273.15);

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
  const auto& air_temp = *S.Get<CompositeVector>(temp_key_, tag).ViewComponent("cell", false);
  const auto& swe = *S.Get<CompositeVector>(snow_key_, tag).ViewComponent("cell", false);
  auto& res = *result[0]->ViewComponent("cell", false);

  if (wrt_key == temp_key_) {
    for (const auto& lc : land_cover_) {
      auto lc_ids = mesh->getSetEntities(
        lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      for (auto c : lc_ids) {
        if (air_temp[0][c] - snow_temp_shift_ > 273.15) {
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
        if (swe[0][c] < lc.second.snow_transition_depth &&
            air_temp[0][c] - snow_temp_shift_ > 273.15) {
          res[0][c] = melt_rate_ * (air_temp[0][c] - snow_temp_shift_ - 273.15) /
                      lc.second.snow_transition_depth;
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
