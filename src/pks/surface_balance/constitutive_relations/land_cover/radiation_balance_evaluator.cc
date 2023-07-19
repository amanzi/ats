/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates a net radiation balance for surface, snow, and canopy.
#include "radiation_balance_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

RadiationBalanceEvaluator::RadiationBalanceEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), compatible_(false)
{
  Key akey = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(akey);
  akey = Keys::getVarName(akey);
  Key dtype = Keys::guessDomainType(domain);
  if (dtype == "surface") {
    domain_surf_ = domain;
    domain_snow_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "snow");
    domain_canopy_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "canopy");
  } else if (dtype == "canopy") {
    domain_canopy_ = domain;
    domain_snow_ = Keys::readDomainHint(plist_, domain_canopy_, "canopy", "snow");
    domain_surf_ = Keys::readDomainHint(plist_, domain_canopy_, "canopy", "surface");
  } else if (dtype == "snow") {
    domain_snow_ = domain;
    domain_canopy_ = Keys::readDomainHint(plist_, domain_snow_, "snow", "canopy");
    domain_surf_ = Keys::readDomainHint(plist_, domain_snow_, "snow", "surface");
  } else {
    domain_surf_ = plist_.get<std::string>("surface domain name");
    domain_snow_ = plist_.get<std::string>("snow domain name");
    domain_canopy_ = plist_.get<std::string>("canopy domain name");
  }
  my_keys_.clear();

  rad_bal_surf_key_ =
    Keys::readKey(plist_, domain_surf_, "surface radiation balance", "radiation_balance");
  my_keys_.emplace_back(KeyTag{ rad_bal_surf_key_, tag });

  rad_bal_snow_key_ =
    Keys::readKey(plist_, domain_snow_, "snow radiation balance", "radiation_balance");
  my_keys_.emplace_back(KeyTag{ rad_bal_snow_key_, tag });

  rad_bal_can_key_ =
    Keys::readKey(plist_, domain_canopy_, "canopy radiation balance", "radiation_balance");
  my_keys_.emplace_back(KeyTag{ rad_bal_can_key_, tag });

  albedo_surf_key_ = Keys::readKey(plist_, domain_surf_, "surface albedos", "albedos");
  dependencies_.insert(KeyTag{ albedo_surf_key_, tag });
  emissivity_surf_key_ =
    Keys::readKey(plist_, domain_surf_, "surface emissivities", "emissivities");
  dependencies_.insert(KeyTag{ emissivity_surf_key_, tag });
  sw_in_key_ = Keys::readKey(
    plist_, domain_surf_, "incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(KeyTag{ sw_in_key_, tag });
  lw_in_key_ = Keys::readKey(
    plist_, domain_surf_, "incoming longwave radiation", "incoming_longwave_radiation");
  dependencies_.insert(KeyTag{ lw_in_key_, tag });

  temp_surf_key_ = Keys::readKey(plist_, domain_surf_, "surface temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_surf_key_, tag });
  temp_snow_key_ = Keys::readKey(plist_, domain_snow_, "snow temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_snow_key_, tag });
  temp_canopy_key_ = Keys::readKey(plist_, domain_canopy_, "canopy temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_canopy_key_, tag });
  area_frac_key_ = Keys::readKey(plist_, domain_surf_, "area fractions", "area_fractions");
  dependencies_.insert(KeyTag{ area_frac_key_, tag });
  lai_key_ = Keys::readKey(plist_, domain_canopy_, "leaf area index", "leaf_area_index");
  dependencies_.insert(KeyTag{ lai_key_, tag });
}


void
RadiationBalanceEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (!compatible_) {
    land_cover_ = getLandCover(S.ICList().sublist("land cover types"),
                               { "beers_k_lw", "beers_k_sw", "albedo_canopy" });

    for (const auto& dep : dependencies_) {
      // dependencies on same mesh, but some have two
      if (dep.first == albedo_surf_key_ || dep.first == emissivity_surf_key_ ||
          dep.first == area_frac_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
          .SetMesh(S.GetMesh(domain_surf_))
          ->SetGhosted(false)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 2);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
          .SetMesh(S.GetMesh(domain_surf_))
          ->SetGhosted(false)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
    }
    compatible_ = true;
  }
}


void
RadiationBalanceEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& rad_bal_surf = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& rad_bal_snow = *results[1]->ViewComponent("cell", false);
  Epetra_MultiVector& rad_bal_can = *results[2]->ViewComponent("cell", false);

  const Epetra_MultiVector& albedo =
    *S.Get<CompositeVector>(albedo_surf_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& emiss =
    *S.Get<CompositeVector>(emissivity_surf_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& sw_in =
    *S.Get<CompositeVector>(sw_in_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& lw_in =
    *S.Get<CompositeVector>(lw_in_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& temp_surf =
    *S.Get<CompositeVector>(temp_surf_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& temp_snow =
    *S.Get<CompositeVector>(temp_snow_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& temp_canopy =
    *S.Get<CompositeVector>(temp_canopy_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& area_frac =
    *S.Get<CompositeVector>(area_frac_key_, tag).ViewComponent("cell", false);
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

      // lw out of each layer
      double lw_surf = Relations::OutgoingLongwaveRadiation(temp_surf[0][c], emiss[0][c]);
      double lw_snow = Relations::OutgoingLongwaveRadiation(temp_snow[0][c], emiss[1][c]);
      double lw_can = Relations::OutgoingLongwaveRadiation(temp_canopy[0][c], e_can_lw);

      // surface connections
      double lw_down = lw_atm_surf + lw_can;
      double lw_up_surf = (1 - emiss[0][c]) * lw_down + lw_surf;
      double lw_up_snow = (1 - emiss[1][c]) * lw_down + lw_snow;

      // radiation balances -- see Figure 4.1 in CLM Tech Note
      rad_bal_surf[0][c] = (1 - albedo[0][c]) * sw_atm_surf + lw_down - lw_up_surf;
      rad_bal_snow[0][c] = (1 - albedo[1][c]) * sw_atm_surf + lw_down - lw_up_snow;

      rad_bal_can[0][c] = (1 - lc.second.albedo_canopy) * sw_atm_can + lw_atm_can - 2 * lw_can +
                          area_frac[0][c] * e_can_lw * lw_up_surf +
                          area_frac[1][c] * e_can_lw * lw_up_snow;
    }
  }
}

void
RadiationBalanceEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& results)
{
  for (const auto& res : results) res->PutScalar(0.);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
