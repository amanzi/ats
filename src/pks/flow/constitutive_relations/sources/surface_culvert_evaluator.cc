/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)
*/

#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "surface_culvert_evaluator.hh"
#include "FunctionFactory.hh"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

// Helper function to compute area-weighted average of a variable over a set of mesh entities
inline double
computeAreaWeightedAverage(
  const AmanziMesh::Mesh& mesh,
  const Amanzi::AmanziMesh::MeshCache<Amanzi::MemSpace_kind::HOST>::cEntity_ID_View& entity_list,
  const Epetra_MultiVector& cv,
  const Epetra_MultiVector& var)
{
  double terms_local[2] = { 0, 0 };
  for (auto c : entity_list) {
    terms_local[0] += cv[0][c] * var[0][c];
    terms_local[1] += cv[0][c];
  }
  double terms_global[2] = { 0, 0 };
  mesh.getComm()->SumAll(terms_local, terms_global, 2);
  return terms_global[0] / terms_global[1]; // Area-weighted average
}


// Constructor: Set up dependencies and culvert parameters
SurfCulvertEvaluator::SurfCulvertEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  // Set up dependencies for required field variables
  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  pd_key_ = Keys::readKey(plist, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });

  wc_key_ = Keys::readKey(plist, domain_, "water content", "water_content");
  dependencies_.insert(KeyTag{ wc_key_, tag });

  pe_key_ = Keys::readKey(plist, domain_, "potential", "pres_elev");
  dependencies_.insert(KeyTag{ pe_key_, tag });

  elev_key_ = Keys::readKey(plist, domain_, "elevation", "elevation");
  dependencies_.insert(KeyTag{ elev_key_, tag });

  liq_den_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ liq_den_key_, tag });

  // Read culvert configuration parameters
  culvert_inlet_region_ = plist.get<std::string>("culvert inlet region");
  culvert_outlet_region_ = plist.get<std::string>("culvert outlet region");
  Nb_ = plist.get<int>("number of barrels", 1);
  L_feet_ = plist.get<double>("culvert length", 10) * 3.28084;
  D_feet_ = plist.get<double>("culvert diameter", 1) * 3.28084;
  n_ = plist.get<double>("culvert roughness coefficient", 0.013);
  C_ = plist.get<double>("culvert discharge coefficient", 0.6);
}

void
SurfCulvertEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  double dt = S.Get<double>("dt", tag);

  // Get required field variables
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& pd = *S.Get<CompositeVector>(pd_key_, tag).ViewComponent("cell", false);
  const auto& liq_den = *S.Get<CompositeVector>(liq_den_key_, tag).ViewComponent("cell", false);
  const auto& wc = *S.Get<CompositeVector>(wc_key_, tag).ViewComponent("cell", false);
  const auto& pe = *S.Get<CompositeVector>(pe_key_, tag).ViewComponent("cell", false);
  const auto& elev = *S.Get<CompositeVector>(elev_key_, tag).ViewComponent("cell", false);

  auto& surf_src = *result[0]->ViewComponent("cell");
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  // Get mesh entities for inlet and outlet regions
  auto inlet_id_list = mesh.getSetEntities(
    culvert_inlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  auto outlet_id_list = mesh.getSetEntities(
    culvert_outlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // Compute area-weighted averages for hydraulic calculations
  double avg_pe_inlet = computeAreaWeightedAverage(mesh, inlet_id_list, cv, pe);
  double avg_elev_inlet = computeAreaWeightedAverage(mesh, inlet_id_list, cv, elev);
  double avg_pe_outlet = computeAreaWeightedAverage(mesh, outlet_id_list, cv, pe);

  // Convert to feet for hydraulic calculations (using US customary units)
  double h_up_feet = avg_pe_inlet * 3.28084;
  double h_down_feet = avg_pe_outlet * 3.28084;
  double z_invert_feet = avg_elev_inlet * 3.28084;

  // Calculate culvert geometry parameters
  double A = M_PI * D_feet_ * D_feet_ / 4.0; // Cross-sectional area
  double P = M_PI * D_feet_;                 // Wetted perimeter
  double R = A / P;                          // Hydraulic radius
  double d = h_up_feet - z_invert_feet;      // Depth at inlet

  // Calculate inlet head for discharge calculations based on different hydraulic regimes
  double d_eps = 0.01; // feet - small depth threshold
  double h_i = 0.0;    // inlet head for discharge calculation
  if (d <= 0) {
    // No water depth - no flow
    h_i = 0.0;
  } else if (d < d_eps) {
    // Very shallow flow - use cubic ramp to avoid numerical issues
    h_i = 0.5 * d * (d / d_eps);
  } else if (d < D_feet_) {
    // Partial flow regime - culvert not full, inlet control
    h_i = 0.5 * d;
  } else {
    // Full flow regime - culvert running full, pressure flow
    h_i = h_up_feet - (z_invert_feet + D_feet_ / 2.0);
  }


  // Calculate flow rates using culvert hydraulics equations
  double g = 9.81;
  double Q_inlet = 0.0;
  if (h_i > 0.0) Q_inlet = Nb_ * C_ * A * std::sqrt(2.0 * g * h_i);

  double h_o = std::max(h_up_feet - h_down_feet, 0.0001);
  double k = 1.5 + (29.0 * n_ * n_ * L_feet_) / std::pow(R, 1.33);
  double Q_outlet = Nb_ * C_ * A * std::sqrt((2.0 * g * h_o) / k);
  double Q_cfs = (Q_inlet * Q_outlet) / std::sqrt(Q_inlet * Q_inlet + Q_outlet * Q_outlet + 1e-12);


  double Q = Q_cfs * 0.028316847; // Convert from ft³/s to m³/s

  // Check water availability constraints
  double available_water_mol = 0.0;
  for (auto c : inlet_id_list) {
    available_water_mol += wc[0][c] * liq_den[0][c] * cv[0][c];
  }
  double total_available_water_mol = 0.0;
  mesh.getComm()->SumAll(&available_water_mol, &total_available_water_mol, 1);

  // Early return if no flow or no water available
  if (Q == 0.0 || total_available_water_mol <= 0.0) {
    for (auto c : inlet_id_list) surf_src[0][c] = 0.0;
    for (auto c : outlet_id_list) surf_src[0][c] = 0.0;
    return;
  }

  // Limit flow rate by available water
  double mols_required = Q * dt;
  if (mols_required > total_available_water_mol && mols_required > 0.0) {
    Q *= total_available_water_mol / mols_required;
  }

  // Distribute source terms at inlet (water removal)
  double sum_wc_l = 0;
  for (auto c : inlet_id_list) sum_wc_l += wc[0][c];
  double sum_wc_g = 0;
  mesh.getComm()->SumAll(&sum_wc_l, &sum_wc_g, 1);
  double safe_sum_wc_g = std::max(sum_wc_g, 1e-8);

  for (auto c : inlet_id_list) {
    double safe_cv = std::max(cv[0][c], 1e-8);
    surf_src[0][c] = -Q * liq_den[0][c] * wc[0][c] / (safe_sum_wc_g * safe_cv);
  }

  // Distribute source terms at outlet (water addition)
  double sum_cv_l = 0;
  for (auto c : outlet_id_list) sum_cv_l += cv[0][c];
  double sum_cv_g = 0;
  mesh.getComm()->SumAll(&sum_cv_l, &sum_cv_g, 1);
  double safe_sum_cv_g = std::max(sum_cv_g, 1e-8);

  for (auto c : outlet_id_list) {
    surf_src[0][c] = Q * liq_den[0][c] / safe_sum_cv_g;
  }
}

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
