/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates a net radiation balance for surface, snow, and canopy.

#include "canopy_radiation_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

CanopyRadiationEvaluator::CanopyRadiationEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), compatible_(false)
{
  Key akey = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  domain_canopy_ = Keys::getDomain(akey);

  // process my keys
  akey = Keys::getVarName(akey);
  my_keys_.clear();

  if (Keys::in("shortwave", akey)) {
    can_down_sw_key_ =
      Keys::readKey(plist_, domain_canopy_, "canopy downward shortwave radiation", akey);
  } else {
    can_down_sw_key_ = Keys::readKey(plist_,
                                     domain_canopy_,
                                     "canopy downward shortwave radiation",
                                     "downward_shortwave_radiation");
  }
  my_keys_.emplace_back(KeyTag{ can_down_sw_key_, tag });

  if (Keys::in("longwave", akey)) {
    can_down_lw_key_ =
      Keys::readKey(plist_, domain_canopy_, "canopy downward longwave radiation", akey);
  } else {
    can_down_lw_key_ = Keys::readKey(
      plist_, domain_canopy_, "canopy downward longwave radiation", "downward_longwave_radiation");
  }
  my_keys_.emplace_back(KeyTag{ can_down_lw_key_, tag });

  if (Keys::in("balance", akey)) {
    rad_bal_can_key_ = Keys::readKey(plist_, domain_canopy_, "canopy radiation balance", akey);
  } else {
    rad_bal_can_key_ =
      Keys::readKey(plist_, domain_canopy_, "canopy radiation balance", "radiation_balance");
  }
  my_keys_.emplace_back(KeyTag{ rad_bal_can_key_, tag });

  // process dependencies
  sw_in_key_ = Keys::readKey(
    plist_, domain_canopy_, "incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(KeyTag{ sw_in_key_, tag });
  lw_in_key_ = Keys::readKey(
    plist_, domain_canopy_, "incoming longwave radiation", "incoming_longwave_radiation");
  dependencies_.insert(KeyTag{ lw_in_key_, tag });

  temp_canopy_key_ = Keys::readKey(plist_, domain_canopy_, "canopy temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_canopy_key_, tag });

  lai_key_ = Keys::readKey(plist_, domain_canopy_, "leaf area index", "leaf_area_index");
  dependencies_.insert(KeyTag{ lai_key_, tag });
}


void
CanopyRadiationEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (!compatible_) {
    land_cover_ =
      getLandCover(S.ICList().sublist("land cover types"),
                   { "beers_k_lw", "beers_k_sw", "emissivity_canopy", "albedo_canopy" });

    for (const auto& dep : dependencies_) {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
        .SetMesh(S.GetMesh(domain_canopy_))
        ->SetGhosted(false)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
    compatible_ = true;
  }
}


void
CanopyRadiationEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& down_sw = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& down_lw = *results[1]->ViewComponent("cell", false);
  Epetra_MultiVector& rad_bal_can = *results[2]->ViewComponent("cell", false);

  const Epetra_MultiVector& sw_in =
    *S.Get<CompositeVector>(sw_in_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& lw_in =
    *S.Get<CompositeVector>(lw_in_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& temp_canopy =
    *S.Get<CompositeVector>(temp_canopy_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& lai =
    *S.Get<CompositeVector>(lai_key_, tag).ViewComponent("cell", false);

  auto mesh = results[0]->Mesh();

  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh->getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : lc_ids) {
      // NOTE: emissivity = absorptivity, we use e to notate both
      // Beer's law to find absorptivity of canopy
      double e_can_sw = Relations::BeersLawAbsorptivity(lc.second.beers_k_sw, lai[0][c]);
      double e_can_lw = Relations::BeersLawAbsorptivity(lc.second.beers_k_lw, lai[0][c]);

      // sw atm to canopy and surface
      double sw_atm_can = e_can_sw * sw_in[0][c];
      double sw_atm_surf = sw_in[0][c] - sw_atm_can;

      // lw atm to canopy and surface
      double lw_atm_can = e_can_lw * lw_in[0][c];
      double lw_atm_surf = lw_in[0][c] - lw_atm_can;

      // black-body radiation for LW out
      double lw_can = Relations::OutgoingLongwaveRadiation(temp_canopy[0][c], e_can_lw);
      down_sw[0][c] = sw_atm_surf;
      down_lw[0][c] = lw_atm_surf + lw_can;

      // factor of 2 is for upward and downward longwave emission NOTE: this is
      // missing surface-to-canopy LW radiation, which is included externally
      // with another evaluator as it is calculated later in the surface
      // balance evaluator.
      rad_bal_can[0][c] = (1 - lc.second.albedo_canopy) * sw_atm_can + lw_atm_can - 2 * lw_can;
    }
  }
}

void
CanopyRadiationEvaluator::EvaluatePartialDerivative_(const State& S,
                                                     const Key& wrt_key,
                                                     const Tag& wrt_tag,
                                                     const std::vector<CompositeVector*>& results)
{
  for (const auto& res : results) res->PutScalar(0.);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
