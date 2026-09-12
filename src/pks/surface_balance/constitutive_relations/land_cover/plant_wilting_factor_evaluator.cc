/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Plant wilting factor provides a moisture availability-based limiter on transpiration.
#include "plant_wilting_factor_evaluator.hh"
#include "plant_wilting_factor_model.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

const std::string PlantWiltingFactorEvaluator::eval_type = "plant wilting factor";

// Constructor from ParameterList
PlantWiltingFactorEvaluator::PlantWiltingFactorEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_keys_.front().first);
  domain_surf_ = Keys::readDomainHint(*plist_, domain_sub_, "domain", "surface");

  // - pull Keys from plist
  pc_key_ = Keys::readKey(*plist_, domain_sub_, "capillary pressure", "capillary_pressure_gas_liq");
  dependencies_.insert(KeyTag{ pc_key_, tag });

  land_cover_ =
    getLandCoverMap(plist_->sublist("model parameters"),
                    { "stomata_closed_capillary_pressure", "stomata_open_capillary_pressure" });
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
PlantWiltingFactorEvaluator::Clone() const
{
  return Teuchos::rcp(new PlantWiltingFactorEvaluator(*this));
}


void
PlantWiltingFactorEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  auto pc_v = S.Get<CompositeVector>(pc_key_, tag).viewComponent("cell", false);
  auto result_v = result[0]->viewComponent("cell", false);

  auto surf_mesh = S.GetMesh(domain_surf_);
  const AmanziMesh::MeshCache& subsurf_mesh = S.GetMesh(domain_sub_)->getCache();

  for (const auto& region_lc : land_cover_) {
    auto lc_ids = surf_mesh->getSetEntities<MemSpace_kind::DEVICE>(
      region_lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    PlantWiltingFactorModel model(region_lc.second);

    Kokkos::parallel_for(
      "PlantWiltingFactorEvaluator::Evaluate", lc_ids.size(), KOKKOS_LAMBDA(const int i) {
        AmanziMesh::Entity_ID sc = lc_ids(i);
        for (auto c : subsurf_mesh.columns.getCells<MemSpace_kind::DEVICE>(sc)) {
          result_v(c, 0) = model.PlantWiltingFactor(pc_v(c, 0));
        }
      });
  }
}


void
PlantWiltingFactorEvaluator::EvaluatePartialDerivative_(const State& S,
                                                        const Key& wrt_key,
                                                        const Tag& wrt_tag,
                                                        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  if (wrt_key == pc_key_) {
    auto pc_v = S.Get<CompositeVector>(pc_key_, tag).viewComponent("cell", false);
    auto result_v = result[0]->viewComponent("cell", false);

    auto surf_mesh = S.GetMesh(domain_surf_);
    const AmanziMesh::MeshCache& subsurf_mesh = S.GetMesh(domain_sub_)->getCache();

    for (const auto& region_lc : land_cover_) {
      auto lc_ids = surf_mesh->getSetEntities<MemSpace_kind::DEVICE>(
        region_lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      PlantWiltingFactorModel model(region_lc.second);

      Kokkos::parallel_for(
        "PlantWiltingFactorEvaluator::EvaluatePartialDerivative",
        lc_ids.size(),
        KOKKOS_LAMBDA(const int i) {
          AmanziMesh::Entity_ID sc = lc_ids(i);
          for (auto c : subsurf_mesh.columns.getCells<MemSpace_kind::DEVICE>(sc)) {
            result_v(c, 0) = model.DPlantWiltingFactorDCapillaryPressureGasLiq(pc_v(c, 0));
          }
        });
    }
  }
}


void
PlantWiltingFactorEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S);
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
