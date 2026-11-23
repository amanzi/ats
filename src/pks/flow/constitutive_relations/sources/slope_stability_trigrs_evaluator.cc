/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

//! Provides slope stability as a factor of safety for landslide.
#include "slope_stability_trigrs_evaluator.hh"
#include "slope_stability_trigrs_model.hh"

namespace Amanzi {
namespace Flow {

// Constructor from ParameterList
SlopeStabilityTrigrsEvaluator::SlopeStabilityTrigrsEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  InitializeFromPlist_();

  std::string params_name = plist_.get<std::string>("model parameters", "WRM parameters");
  Teuchos::ParameterList& sublist = plist_.sublist(params_name);
  models_ = createRegionModelPartition<SlopeStabilityTrigrsModel>(sublist);
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
SlopeStabilityTrigrsEvaluator::Clone() const
{
  return Teuchos::rcp(new SlopeStabilityTrigrsEvaluator(*this));
}


// Initialize by setting up dependencies
void
SlopeStabilityTrigrsEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_keys_.front().first);
  domain_surf_ = Keys::readDomainHint(plist_, domain_sub_, "domain", "surface");
  Tag tag = my_keys_.front().second;

  internal_friction_key_ = Keys::readKey(plist_, domain_sub_,
    "effective internal friction", "effective_internal_friction");
  dependencies_.insert(KeyTag{ internal_friction_key_, tag });

  soil_cohesion_key_ = Keys::readKey(plist_, domain_sub_,
    "effective soil cohesion", "effective_soil_cohesion");
  dependencies_.insert(KeyTag{ soil_cohesion_key_, tag });

  pressure_key_ = Keys::readKey(plist_, domain_sub_, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pressure_key_, tag });

  sat_key_ = Keys::readKey(plist_, domain_sub_, "saturation", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_key_, tag });

  depth_key_ = Keys::readKey(plist_, domain_sub_, "depth", "depth");
  dependencies_.insert(KeyTag{ depth_key_, tag });

  rho_r_key_ = Keys::readKey(plist_, domain_sub_, "density rock", "density_rock");
  dependencies_.insert(KeyTag{ rho_r_key_, tag });

  poro_key_ = Keys::readKey(plist_, domain_sub_, "porosity", "porosity");
  dependencies_.insert(KeyTag{ poro_key_, tag });

  rhow_key_ = Keys::readKey(plist_, domain_sub_, "mass density", "mass_density_liquid");
  dependencies_.insert(KeyTag{ rhow_key_, tag });

  slope_key_ = Keys::readKey(plist_, domain_surf_, "slope", "slope_magnitude");
  dependencies_.insert(KeyTag{ slope_key_, tag });
}


void
SlopeStabilityTrigrsEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result[0]->Mesh(), -1);
    models_->first->Verify();
  }

  const Epetra_MultiVector& internal_friction_v =
    *S.Get<CompositeVector>(internal_friction_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& soil_cohesion_v =
    *S.Get<CompositeVector>(soil_cohesion_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& pressure_v =
    *S.Get<CompositeVector>(pressure_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& sat_v =
    *S.Get<CompositeVector>(sat_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& depth_v =
    *S.Get<CompositeVector>(depth_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& rho_r_v =
    *S.Get<CompositeVector>(rho_r_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& poro_v =
    *S.Get<CompositeVector>(poro_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& rhow_v =
    *S.Get<CompositeVector>(rhow_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& slope_v =
    *S.Get<CompositeVector>(slope_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  double p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  double gravity = -(S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT))[2];

  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
    double sigma = std::atan(slope_v[0][col]); // convert to rad
    const auto& col_cells = mesh.columns.getCells(col);
    for (int i = 0; i != col_cells.size(); ++i) {
      AmanziMesh::Entity_ID c = col_cells[i];
      int region = (*models_->first)[c];
      const auto& model = *models_->second[region];
      result_v[0][c] = model.SlopeStabilityByTrigrs(
        internal_friction_v[0][c], sigma, soil_cohesion_v[0][c],
        pressure_v[0][c] - p_atm, sat_v[0][c], rho_r_v[0][c],
        depth_v[0][c], gravity, poro_v[0][c], rhow_v[0][c]);
    }
  }
}


void
SlopeStabilityTrigrsEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result[0]->Mesh(), -1);
    models_->first->Verify();
  }

  const Epetra_MultiVector& internal_friction_v =
    *S.Get<CompositeVector>(internal_friction_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& soil_cohesion_v =
    *S.Get<CompositeVector>(soil_cohesion_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& pressure_v =
    *S.Get<CompositeVector>(pressure_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& sat_v =
    *S.Get<CompositeVector>(sat_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& depth_v =
    *S.Get<CompositeVector>(depth_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& rho_r_v =
    *S.Get<CompositeVector>(rho_r_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& poro_v =
    *S.Get<CompositeVector>(poro_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& rhow_v =
    *S.Get<CompositeVector>(rhow_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& slope_v =
    *S.Get<CompositeVector>(slope_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  double p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);
  double gravity = -(S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT))[2];

  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  if (wrt_key == pressure_key_) {
    for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
      double sigma = std::atan(slope_v[0][col]); // convert to rad
      const auto& col_cells = mesh.columns.getCells(col);
      for (int i = 0; i != col_cells.size(); ++i) {
        AmanziMesh::Entity_ID c = col_cells[i];
        int region = (*models_->first)[c];
        const auto& model = *models_->second[region];
        result_v[0][c] = model.DSlopeStabilityByTrigrsDpressure(
          internal_friction_v[0][c], sigma, sat_v[0][c],
          rho_r_v[0][c], depth_v[0][c], gravity, 
          poro_v[0][c], rhow_v[0][c]);
      }
    }
  } else if (wrt_key == sat_key_) {
    for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
      double sigma = std::atan(slope_v[0][col]); // convert to rad
      const auto& col_cells = mesh.columns.getCells(col);
      for (int i = 0; i != col_cells.size(); ++i) {
        AmanziMesh::Entity_ID c = col_cells[i];
        int region = (*models_->first)[c];
        const auto& model = *models_->second[region];
        result_v[0][c] = model.DSlopeStabilityByTrigrsDSat(
          internal_friction_v[0][c], sigma, soil_cohesion_v[0][c],
          pressure_v[0][c] - p_atm, sat_v[0][c], rho_r_v[0][c],
          depth_v[0][c], gravity, poro_v[0][c], rhow_v[0][c]);
      }
    }
  } else if (wrt_key == poro_key_) {
    for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
      double sigma = std::atan(slope_v[0][col]); // convert to rad
      const auto& col_cells = mesh.columns.getCells(col);
      for (int i = 0; i != col_cells.size(); ++i) {
        AmanziMesh::Entity_ID c = col_cells[i];
        int region = (*models_->first)[c];
        const auto& model = *models_->second[region];
        result_v[0][c] = model.DSlopeStabilityByTrigrsDPorosity(
          internal_friction_v[0][c], sigma, soil_cohesion_v[0][c],
          pressure_v[0][c] - p_atm, sat_v[0][c], rho_r_v[0][c],
          depth_v[0][c], gravity, poro_v[0][c], rhow_v[0][c]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}


void
SlopeStabilityTrigrsEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // dependencies on the subsurface mesh
  CompositeVectorSpace dep_fac;
  dep_fac.SetMesh(S.GetMesh(domain_sub_))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // dependencies on the surface mesh
  CompositeVectorSpace dep_fac_surf;
  dep_fac_surf.SetMesh(S.GetMesh(domain_surf_))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  for (const auto& dep : dependencies_) {
    if (Keys::getDomain(dep.first) == domain_sub_) {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(dep_fac);
    } else {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(dep_fac_surf);
    }
  }
}

} // namespace Flow
} // namespace Amanzi
