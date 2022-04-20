/*
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
  : SecondaryVariablesFieldEvaluator(plist),
    compatible_(false)
{
  Key domain = Keys::getDomain(Keys::cleanPListName(plist_.name()));
  Key dtype = Keys::guessDomainType(domain);
  if (dtype == "surface") {
    domain_surf_ = domain;
    domain_snow_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "snow");
    domain_canopy_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "canopy");
  } else if (dtype == "canopy") {
    domain_canopy_ = domain;
    domain_snow_ = Keys::readDomainHint(plist_, domain_canopy_, "canopy", "snow");
    domain_surf_ = Keys::readDomainHint(plist_, domain_canopy_, "canopy", "surface");
  } else {
    domain_surf_ = plist_.get<std::string>("surface domain name");
    domain_snow_ = plist_.get<std::string>("snow domain name");
    domain_canopy_ = plist_.get<std::string>("canopy domain name");
  }

  rad_bal_surf_key_ = Keys::readKey(plist_, domain_surf_, "surface radiation balance", "radiation_balance");
  rad_bal_snow_key_ = Keys::readKey(plist_, domain_snow_, "snow radiation balance", "radiation_balance");
  rad_bal_can_key_ = Keys::readKey(plist_, domain_canopy_, "canopy radiation balance", "radiation_balance");

  my_keys_.push_back(rad_bal_surf_key_);
  my_keys_.push_back(rad_bal_snow_key_);
  my_keys_.push_back(rad_bal_can_key_);

  albedo_surf_key_ = Keys::readKey(plist_, domain_surf_, "surface albedos", "albedos");
  dependencies_.insert(albedo_surf_key_);
  emissivity_surf_key_ = Keys::readKey(plist_, domain_surf_, "surface emissivities", "emissivities");
  dependencies_.insert(emissivity_surf_key_);
  sw_in_key_ = Keys::readKey(plist_, domain_surf_, "incoming shortwave radiation", "incoming_shortwave_radiation");
  dependencies_.insert(sw_in_key_);
  lw_in_key_ = Keys::readKey(plist_, domain_surf_, "incoming longwave radiation", "incoming_longwave_radiation");
  dependencies_.insert(lw_in_key_);

  temp_surf_key_ = Keys::readKey(plist_, domain_surf_, "surface temperature", "temperature");
  dependencies_.insert(temp_surf_key_);
  temp_snow_key_ = Keys::readKey(plist_, domain_snow_, "snow temperature", "temperature");
  dependencies_.insert(temp_snow_key_);
  temp_canopy_key_ = Keys::readKey(plist_, domain_canopy_, "canopy temperature", "temperature");
  dependencies_.insert(temp_canopy_key_);
  area_frac_key_ = Keys::readKey(plist_, domain_surf_, "area fractions", "area_fractions");
  dependencies_.insert(area_frac_key_);
  lai_key_ = Keys::readKey(plist_, domain_canopy_, "leaf area index", "leaf_area_index");
  dependencies_.insert(lai_key_);
}

void
RadiationBalanceEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (!compatible_) {
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
            {"beers_law_lw", "beers_law_sw", "emissivity_canopy", "albedo_canopy"});

    for (const auto& my_key : my_keys_) {
      // require all domains are the same
      S->RequireField(my_key, my_key)
        ->SetMesh(S->GetMesh(domain_surf_))
        ->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }

    for (const auto& dep : dependencies_) {
      // dependencies on same mesh, but some have two
      if (dep == albedo_surf_key_ ||
          dep == emissivity_surf_key_ ||
          dep == area_frac_key_) {
        S->RequireField(dep)
          ->SetMesh(S->GetMesh(domain_surf_))
          ->SetGhosted(false)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 2);
      } else {
        S->RequireField(dep)
          ->SetMesh(S->GetMesh(domain_surf_))
          ->SetGhosted(false)
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
      // Recurse into the tree to propagate info to leaves.
      S->RequireFieldEvaluator(dep)->EnsureCompatibility(S);
    }
    compatible_ = true;
  }
}

void
RadiationBalanceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
               const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  Epetra_MultiVector& rad_bal_surf = *results[0]->ViewComponent("cell",false);
  Epetra_MultiVector& rad_bal_snow = *results[1]->ViewComponent("cell",false);
  Epetra_MultiVector& rad_bal_can = *results[2]->ViewComponent("cell",false);

  const Epetra_MultiVector& albedo = *S->GetFieldData(albedo_surf_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& emiss = *S->GetFieldData(emissivity_surf_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& sw_in = *S->GetFieldData(sw_in_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& lw_in = *S->GetFieldData(lw_in_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& temp_surf = *S->GetFieldData(temp_surf_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& temp_snow = *S->GetFieldData(temp_snow_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& temp_canopy = *S->GetFieldData(temp_canopy_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& area_frac = *S->GetFieldData(area_frac_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& lai = *S->GetFieldData(lai_key_)->ViewComponent("cell",false);

  auto mesh = results[0]->Mesh();

  for (const auto& lc : land_cover_) {
    AmanziMesh::Entity_ID_List lc_ids;
    mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    for (auto c : lc_ids) {
      // Beer's law to find attenuation of radiation to surface in sw
      double sw_atm_surf = Relations::BeersLaw(sw_in[0][c], lc.second.beers_k_sw, lai[0][c]);
      double sw_atm_can = sw_in[0][c] - sw_atm_surf;

      // Beer's law to find attenuation of radiation to surface in lw -- note
      // this should be almost 0 for any LAI
      double lw_atm_surf = Relations::BeersLaw(lw_in[0][c], lc.second.beers_k_lw, lai[0][c]);
      double lw_atm_can = lw_in[0][c] - lw_atm_surf;

      // black-body radiation for LW out
      double lw_surf = Relations::OutgoingLongwaveRadiation(temp_surf[0][c], emiss[0][c]);
      double lw_snow = Relations::OutgoingLongwaveRadiation(temp_snow[0][c], emiss[1][c]);
      double lw_can = Relations::OutgoingLongwaveRadiation(temp_canopy[0][c],
              lc.second.emissivity_canopy);

      rad_bal_surf[0][c] = (1-albedo[0][c])*sw_atm_surf + lw_atm_surf + lw_can - lw_surf;
      rad_bal_snow[0][c] = (1-albedo[1][c])*sw_atm_surf + lw_atm_surf + lw_can - lw_snow;
      rad_bal_can[0][c] = (1-lc.second.albedo_canopy)*sw_atm_can + lw_atm_can
        + area_frac[1][c]*lw_snow + area_frac[0][c]*lw_surf - 2*lw_can;
    }
  }
}

void
RadiationBalanceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results)
{
  for (const auto& res : results) res->PutScalar(0.);
//  AMANZI_ASSERT(false);
}


}  // namespace Relations
}  // namespace SurfaceBalance
}  // namespace Amanzi
