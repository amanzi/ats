/*
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

// Constructor from ParameterList
PlantWiltingFactorEvaluator::PlantWiltingFactorEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_key_);
  domain_surf_ = Keys::readDomainHint(plist_, domain_sub_, "domain", "surface");

  // - pull Keys from plist
  pc_key_ = Keys::readKey(plist_, domain_sub_, "capillary pressure", "capillary_pressure_gas_liq");
  dependencies_.insert(pc_key_);
}

// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
PlantWiltingFactorEvaluator::Clone() const
{
  return Teuchos::rcp(new PlantWiltingFactorEvaluator(*this));
}



void
PlantWiltingFactorEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  if (models_.size() == 0) {
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
            {"stomata_closed_mafic_potential", "stomata_open_mafic_potential"});
    for (const auto& lc : land_cover_) {
      models_[lc.first] = Teuchos::rcp(new PlantWiltingFactorModel(lc.second));
    }
  }

  const Epetra_MultiVector& pc_v = *S->GetFieldData(pc_key_)->ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  auto& subsurf_mesh = *S->GetMesh(domain_sub_);
  auto& surf_mesh = *S->GetMesh(domain_surf_);

  for (const auto& region_model : models_) {
    AmanziMesh::Entity_ID_List lc_ids;
    surf_mesh.get_set_entities(region_model.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    for (int sc : lc_ids) {
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        result_v[0][c] = region_model.second->PlantWiltingFactor(pc_v[0][c]);
      }
    }
  }
}


void
PlantWiltingFactorEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  if (wrt_key == pc_key_) {
    const Epetra_MultiVector& pc_v = *S->GetFieldData(pc_key_)->ViewComponent("cell", false);
    Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

    auto& subsurf_mesh = *S->GetMesh(domain_sub_);
    auto& surf_mesh = *S->GetMesh(domain_surf_);

    for (const auto& region_model : models_) {
      AmanziMesh::Entity_ID_List lc_ids;
      surf_mesh.get_set_entities(region_model.first, AmanziMesh::Entity_kind::CELL,
              AmanziMesh::Parallel_type::OWNED, &lc_ids);

      for (int sc : lc_ids) {
        for (auto c : subsurf_mesh.cells_of_column(sc)) {
          result_v[0][c] = region_model.second->DPlantWiltingFactorDCapillaryPressureGasLiq(pc_v[0][c]);
        }
      }
    }
  }
}


} //namespace
} //namespace
} //namespace
