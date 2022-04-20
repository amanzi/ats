/*
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
EvaporationDownregulationEvaluator::EvaporationDownregulationEvaluator(Teuchos::ParameterList& plist) :
  SecondaryVariableFieldEvaluator(plist),
  consistent_(false)
{
  InitializeFromPlist_();
}

// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
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
  domain_surf_ = Keys::getDomain(my_key_);
  domain_sub_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "domain");

  // sat gas and porosity on subsurface
  sat_gas_key_ = Keys::readKey(plist_, domain_sub_, "saturation gas", "saturation_gas");
  dependencies_.insert(sat_gas_key_);

  poro_key_ = Keys::readKey(plist_, domain_sub_, "porosity", "porosity");
  dependencies_.insert(poro_key_);

  // dependency: potential_evaporation on surface
  pot_evap_key_ = Keys::readKey(plist_, domain_surf_, "potential evaporation", "potential_evaporation");
  dependencies_.insert(pot_evap_key_);
}


void
EvaporationDownregulationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& sat_gas = *S->GetFieldData(sat_gas_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S->GetFieldData(poro_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& pot_evap = *S->GetFieldData(pot_evap_key_)->ViewComponent("cell",false);
  Epetra_MultiVector& surf_evap = *result->ViewComponent("cell",false);
  auto& sub_mesh = *S->GetMesh(domain_sub_);
  auto& surf_mesh = *S->GetMesh(domain_surf_);

  for (const auto& region_model : models_) {
    AmanziMesh::Entity_ID_List lc_ids;
    surf_mesh.get_set_entities(region_model.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    for (AmanziMesh::Entity_ID sc : lc_ids) {
      auto c = sub_mesh.cells_of_column(sc)[0];
      surf_evap[0][sc] = region_model.second->Evaporation(sat_gas[0][c], poro[0][c], pot_evap[0][sc]);
    }
  }
}


void
EvaporationDownregulationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  if (wrt_key == pot_evap_key_) {
    const Epetra_MultiVector& sat_gas = *S->GetFieldData(sat_gas_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& poro = *S->GetFieldData(poro_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& pot_evap = *S->GetFieldData(pot_evap_key_)->ViewComponent("cell",false);
    Epetra_MultiVector& surf_evap = *result->ViewComponent("cell",false);
    auto& sub_mesh = *S->GetMesh(domain_sub_);
    auto& surf_mesh = *S->GetMesh(domain_surf_);

    for (const auto& region_model : models_) {
      AmanziMesh::Entity_ID_List lc_ids;
      surf_mesh.get_set_entities(region_model.first, AmanziMesh::Entity_kind::CELL,
              AmanziMesh::Parallel_type::OWNED, &lc_ids);
      for (AmanziMesh::Entity_ID sc : lc_ids) {
        auto c = sub_mesh.cells_of_column(sc)[0];
        surf_evap[0][sc] = region_model.second->DEvaporationDPotentialEvaporation(sat_gas[0][c], poro[0][c], pot_evap[0][sc]);
      }
    }
  }
}


void
EvaporationDownregulationEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (!consistent_) {
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
            {"dessicated_zone_thickness", "clapp_horn_b"});
    for (const auto& lc : land_cover_) {
      models_[lc.first] = Teuchos::rcp(new EvaporationDownregulationModel(lc.second));
    }

    S->RequireField(my_key_, my_key_)
      ->SetMesh(S->GetMesh(domain_surf_))
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    S->RequireField(poro_key_)
      ->SetMesh(S->GetMesh(domain_sub_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(sat_gas_key_)
      ->SetMesh(S->GetMesh(domain_sub_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireField(pot_evap_key_)
      ->SetMesh(S->GetMesh(domain_surf_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    // Check plist for vis or checkpointing control.
    bool io_my_key = plist_.get<bool>("visualize", true);
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

    // Recurse into the tree to propagate info to leaves.
    // Loop over dependencies, making sure they are the same mesh
    for (auto key : dependencies_) {
      S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
    }

    consistent_ = true;
  }
}


} //namespace
} //namespace
} //namespace
