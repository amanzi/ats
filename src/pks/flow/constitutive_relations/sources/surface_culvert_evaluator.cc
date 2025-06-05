/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: 
*/

#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "surface_culvert_evaluator.hh"
#include "FunctionFactory.hh"
#include <cmath>
#include <algorithm>

namespace Amanzi {
namespace Flow {
namespace Relations {

// Helper function to compute area-weighted average
inline double computeAreaWeightedAverage(const AmanziMesh::Mesh& mesh,
                                  const Amanzi::AmanziMesh::MeshCache<Amanzi::MemSpace_kind::HOST>::cEntity_ID_View& entity_list,
                                  const Epetra_MultiVector& cv,
                                  const Epetra_MultiVector& var) {
    double terms_local[2] = { 0, 0 };
    for (auto c : entity_list) {
        terms_local[0] += cv[0][c] * var[0][c];
        terms_local[1] += cv[0][c];
    }
    double terms_global[2] = { 0, 0 };
    mesh.getComm()->SumAll(terms_local, terms_global, 2);
    return terms_global[0] / terms_global[1]; // Area-weighted average
}


SurfCulvertEvaluator::SurfCulvertEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  // dependencies
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

  culvert_inlet_region_ = plist.get<std::string>("culvert inlet region");
  culvert_outlet_region_ = plist.get<std::string>("culvert outlet region");
  //culvert geometry
  Nb_ = plist.get<int>("number of barrels", 1);
  L_ = plist.get<double>("culvert length", 10);
  D_ = plist.get<double>("culvert diameter", 1);
  n_ = plist.get<double>("culvert roughness coefficient", 0.013);
  C_ = plist.get<double>("culvert discharge coefficient", 0.6);

}

// Required methods from SecondaryVariableFieldEvaluator

void
SurfCulvertEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result) 
{

  Tag tag = my_keys_.front().second;
  
  double dt = S.Get<double>("dt", tag);
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& pd = *S.Get<CompositeVector>(pd_key_, tag).ViewComponent("cell", false);
  const auto& liq_den = *S.Get<CompositeVector>(liq_den_key_, tag).ViewComponent("cell", false);
  const auto& wc = *S.Get<CompositeVector>(wc_key_, tag).ViewComponent("cell", false);
  const auto& pe = *S.Get<CompositeVector>(pe_key_, tag).ViewComponent("cell", false);
  const auto& elev = *S.Get<CompositeVector>(elev_key_, tag).ViewComponent("cell", false);

  auto& surf_src= *result[0]->ViewComponent("cell"); // not being reference

  double total = 0.0;
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  auto inlet_id_list = mesh.getSetEntities(
    culvert_inlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  auto outlet_id_list = mesh.getSetEntities(
    culvert_outlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // Calculate the flow rate from culvert (from inlet to outlet)
  double Q;
  double avg_pe_inlet = computeAreaWeightedAverage(mesh, inlet_id_list, cv, pe);
  double avg_elev_inlet = computeAreaWeightedAverage(mesh, inlet_id_list, cv, elev);
  double avg_pe_outlet = computeAreaWeightedAverage(mesh, outlet_id_list, cv, pe);

  // Cross-section
  double A = M_PI * D_ * D_ / 4.0;
  double P = M_PI * D_;
  double R = A / P;

  // Flow depth at the invert
  double d = avg_pe_inlet - avg_elev_inlet;

  // Inlet head based on partial or full submergence
  double h_i = 0.0;
  if (d <= 0) {
    h_i = 0.0; // dry
  } else if (d < D_) {
    h_i = 0.5 * d; // partial
  } else {
    // For full submergence, need invert elevation and upstream head
    h_i = avg_pe_inlet - (avg_elev_inlet + D_ / 2.0); // full
  }

  // Inlet control
  double C_inlet = 0.6;
  double g = 9.81;
  double Q_inlet = 0.0;
  if (h_i > 0.0) {
    Q_inlet = Nb_ * C_inlet * A * std::sqrt(2.0 * g * h_i);
  }

  // Outlet control
  double h_o = std::max(avg_pe_inlet - avg_pe_outlet, 0.0001);

  double k = 1.5 + (29.0 * n_ * n_ * L_) /
                    std::pow(R, 1.33);

  // Q_outlet = N_b * A * sqrt((2 * g * h_o) / k)

  // Use C_ (from plist) as the discharge coefficient.
  double Q_outlet = Nb_ * C_ * A * std::sqrt((2.0 * g * h_o) / k);

  // taking the flow with maximum constraint
  Q = std::min(Q_inlet, Q_outlet);
  // if stability is an issue, we can use the following:
  // double Q = (Q_inlet * Q_outlet) / std::sqrt(Q_inlet * Q_inlet + Q_outlet * Q_outlet + 1e-12);

  // Sink to the inlet cells
  double sum_wc_l = 0;
  for (auto c : inlet_id_list) {
    sum_wc_l += wc[0][c];
  }
  double sum_wc_g = 0;
  mesh.getComm()->SumAll(&sum_wc_l, &sum_wc_g, 1);

  for (auto c : inlet_id_list) {
    if (sum_wc_g > 0) {
      surf_src[0][c] = - Q * liq_den[0][c] * wc[0][c]/(sum_wc_g * cv[0][c]); // mol/(m^2 * s)
    }
  }

  // Source to the outlet cells
  double sum_cv_l = 0; 
  for (auto c : outlet_id_list) {
    sum_cv_l += cv[0][c];
  }
  double sum_cv_g = 0;
  mesh.getComm()->SumAll(&sum_cv_l, &sum_cv_g, 1);

  for (auto c : outlet_id_list) {
    surf_src[0][c] = Q * liq_den[0][c] / (sum_cv_g); // mol/(m^2 * s)
  }

}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi 