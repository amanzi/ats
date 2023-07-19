/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Downregulates bare soil evaporation through a dessicated zone.
#include "evaporation_downregulation_evaluator.hh"
#include "evaporation_downregulation_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
EvaporationDownregulationEvaluator::EvaporationDownregulationEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), consistent_(false)
{
  InitializeFromPlist_();
}

// Virtual copy constructor
Teuchos::RCP<Evaluator>
EvaporationDownregulationEvaluator::Clone() const
{
  return Teuchos::rcp(new EvaporationDownregulationEvaluator(*this));
}


// Initialize by setting up dependencies
void
EvaporationDownregulationEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Tag tag = my_keys_.front().second;
  domain_surf_ = Keys::getDomain(my_keys_.front().first);
  domain_sub_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "domain");

  // sat gas, sat liquid, and porosity on subsurface
  sat_gas_key_ = Keys::readKey(plist_, domain_sub_, "saturation gas", "saturation_gas");
  dependencies_.insert(KeyTag{ sat_gas_key_, tag });
  sat_liq_key_ = Keys::readKey(plist_, domain_sub_, "saturation liquid", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_liq_key_, tag });
  poro_key_ = Keys::readKey(plist_, domain_sub_, "porosity", "porosity");
  dependencies_.insert(KeyTag{ poro_key_, tag });

  // dependency: potential_evaporation on surface
  pot_evap_key_ =
    Keys::readKey(plist_, domain_surf_, "potential evaporation", "potential_evaporation");
  dependencies_.insert(KeyTag{ pot_evap_key_, tag });
}


void
EvaporationDownregulationEvaluator::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& sat_gas =
    *S.Get<CompositeVector>(sat_gas_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& sat_liq =
    *S.Get<CompositeVector>(sat_liq_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& poro =
    *S.Get<CompositeVector>(poro_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& pot_evap =
    *S.Get<CompositeVector>(pot_evap_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& surf_evap = *result[0]->ViewComponent("cell", false);
  auto& sub_mesh = *S.GetMesh(domain_sub_);
  auto& surf_mesh = *S.GetMesh(domain_surf_);

  for (const auto& region_model : models_) {
    auto lc_ids = surf_mesh.getSetEntities(
      region_model.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (AmanziMesh::Entity_ID sc : lc_ids) {
      auto c = sub_mesh.columns.getCells(sc)[0];
      surf_evap[0][sc] =
        region_model.second->Evaporation(sat_gas[0][c], poro[0][c], pot_evap[0][sc], sat_liq[0][c]);
    }
  }
}


void
EvaporationDownregulationEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  if (wrt_key == pot_evap_key_) {
    const Epetra_MultiVector& sat_gas =
      *S.Get<CompositeVector>(sat_gas_key_, tag).ViewComponent("cell", false);
    const Epetra_MultiVector& sat_liq =
      *S.Get<CompositeVector>(sat_liq_key_, tag).ViewComponent("cell", false);
    const Epetra_MultiVector& poro =
      *S.Get<CompositeVector>(poro_key_, tag).ViewComponent("cell", false);
    const Epetra_MultiVector& pot_evap =
      *S.Get<CompositeVector>(pot_evap_key_, tag).ViewComponent("cell", false);
    Epetra_MultiVector& surf_evap = *result[0]->ViewComponent("cell", false);
    auto& sub_mesh = *S.GetMesh(domain_sub_);
    auto& surf_mesh = *S.GetMesh(domain_surf_);

    for (const auto& region_model : models_) {
      auto lc_ids = surf_mesh.getSetEntities(region_model.first,
                                 AmanziMesh::Entity_kind::CELL,
                                 AmanziMesh::Parallel_kind::OWNED);
      for (AmanziMesh::Entity_ID sc : lc_ids) {
        auto c = sub_mesh.columns.getCells(sc)[0];
        surf_evap[0][sc] = region_model.second->DEvaporationDPotentialEvaporation(
          sat_gas[0][c], poro[0][c], pot_evap[0][sc], sat_liq[0][c]);
      }
    }
  }
}


void
EvaporationDownregulationEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (!consistent_) {
    land_cover_ = getLandCover(S.ICList().sublist("land cover types"),
                               { "dessicated_zone_thickness", "clapp_horn_b" });
    for (const auto& lc : land_cover_) {
      models_[lc.first] = Teuchos::rcp(new EvaporationDownregulationModel(lc.second));
    }

    Tag tag = my_keys_.front().second;
    S.Require<CompositeVector, CompositeVectorSpace>(poro_key_, tag)
      .SetMesh(S.GetMesh(domain_sub_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(sat_gas_key_, tag)
      .SetMesh(S.GetMesh(domain_sub_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(sat_liq_key_, tag)
      .SetMesh(S.GetMesh(domain_sub_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.Require<CompositeVector, CompositeVectorSpace>(pot_evap_key_, tag)
      .SetMesh(S.GetMesh(domain_surf_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    consistent_ = true;
  }
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
