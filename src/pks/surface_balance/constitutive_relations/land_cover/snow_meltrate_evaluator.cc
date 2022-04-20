/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ahmad Jan (jana@ornl.gov)
*/

#include "Key.hh"
#include "snow_meltrate_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {


SnowMeltRateEvaluator::SnowMeltRateEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    compatibility_checked_(false)
{
  melt_rate_ = plist.get<double>("snow melt rate [mm day^-1 C^-1]", 2.74) * 0.001 / 86400.; // convert mm/day to m/s
  snow_temp_shift_ = plist.get<double>("air-snow temperature difference [C]", 2.0); // snow is typically a few degrees colder than air at melt time

  domain_ = Keys::getDomain(my_key_);
  domain_surf_ = Keys::readDomainHint(plist_, domain_, "snow", "surface");

  temp_key_ = Keys::readKey(plist, domain_surf_, "air temperature", "air_temperature");
  dependencies_.insert(temp_key_);

  snow_key_ = Keys::readKey(plist, domain_, "snow water equivalent", "water_equivalent");
  dependencies_.insert(snow_key_);
}

void
SnowMeltRateEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // new state!
  land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
                             {"snow_transition_depth"});
  SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
}

// Required methods from SecondaryVariableFieldEvaluator
void
SnowMeltRateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto mesh = S->GetMesh(domain_);

  const auto& air_temp = *S->GetFieldData(temp_key_)->ViewComponent("cell", false);
  const auto& swe = *S->GetFieldData(snow_key_)->ViewComponent("cell", false);
  auto& res = *result->ViewComponent("cell", false);

  for (const auto& lc : land_cover_) {
    AmanziMesh::Entity_ID_List lc_ids;
    mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

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

// Required methods from SecondaryVariableFieldEvaluator
void
SnowMeltRateEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  auto mesh = S->GetMesh(domain_);
  const auto& air_temp = *S->GetFieldData(temp_key_)->ViewComponent("cell", false);
  const auto& swe = *S->GetFieldData(snow_key_)->ViewComponent("cell", false);
  auto& res = *result->ViewComponent("cell", false);

  if (wrt_key == temp_key_) {
    for (const auto& lc : land_cover_) {
      AmanziMesh::Entity_ID_List lc_ids;
      mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                             AmanziMesh::Parallel_type::OWNED, &lc_ids);
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
      AmanziMesh::Entity_ID_List lc_ids;
      mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                             AmanziMesh::Parallel_type::OWNED, &lc_ids);
      for (auto c : lc_ids) {
        if (swe[0][c] < lc.second.snow_transition_depth && air_temp[0][c] - snow_temp_shift_ > 273.15) {
          res[0][c] = melt_rate_ * (air_temp[0][c] - snow_temp_shift_ - 273.15) / lc.second.snow_transition_depth;
        } else {
          res[0][c] = 0.0;
        }
      }
    }
  }
}

} //namespace
} //namespace
} //namespace

