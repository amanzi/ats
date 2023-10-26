
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
#include "surface_pump_system_evaluator.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

//** // Helper function to compute area-weighted average
// double computeAreaWeightedAverage(const AmanziMesh::Mesh& mesh, 
//                                   const std::vector<int>& entity_list,
//                                   const CompositeVector& cv,
//                                   const CompositeVector& var) {
//     double terms_local[2] = { 0, 0 }; 
//     for (auto c : entity_list) {
//         terms_local[0] += cv[0][c] * var[0][c];
//         terms_local[1] += cv[0][c];
//     }
//     double terms_global[2] = { 0, 0 };
//     mesh.getComm()->SumAll(terms_local, terms_global, 2);
//     return terms_global[0] / terms_global[1]; // Area-weighted average
// }
//**


SurfPumpEvaluator::SurfPumpEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  pd_key_ = Keys::readKey(plist, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });
  pe_key_ = Keys::readKey(plist, domain_, "pressure elevation", "pres_elev");
  dependencies_.insert(KeyTag{ pe_key_, tag });
  wc_key_ = Keys::readKey(plist, domain_, "water content", "water_content");
  dependencies_.insert(KeyTag{ wc_key_, tag });
  liq_den_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  pump_outlet_region_ = plist.get<std::string>("pump outlet region", "");
  pump_inlet_region_ = plist.get<std::string>("pump inlet region");
  on_off_region_ = plist.get<std::string>("on off reference region", "");
  max_elev_pumpline_ = plist.get<double>("maximum pumpline elevation", -9999);
  stage_on_ = plist.get<double>("pump start at stage", -9999);
  stage_off_ = plist.get<double>("pump stop at stage", -9999);

  Teuchos::ParameterList& pump_func = plist.sublist("function");
  FunctionFactory fac;
  Q_pump_ = Teuchos::rcp(fac.Create(pump_func));
}

// Required methods from SecondaryVariableFieldEvaluator
void
SurfPumpEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result) 
{

  Tag tag = my_keys_.front().second;
  
  double dt = S.Get<double>("dt", tag);
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& pd = *S.Get<CompositeVector>(pd_key_, tag).ViewComponent("cell", false);
  const auto& liq_den = *S.Get<CompositeVector>(liq_den_key_, tag).ViewComponent("cell", false);
  const auto& wc = *S.Get<CompositeVector>(wc_key_, tag).ViewComponent("cell", false);
  const auto& pe = *S.Get<CompositeVector>(pe_key_, tag).ViewComponent("cell", false);
  

  auto& surf_src= *result[0]->ViewComponent("cell"); // not being reference

  double total = 0.0;
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  auto pump_inlet_id_list = mesh.getSetEntities(
    pump_inlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  std::vector<int> pump_outlet_id_list;
  if (!pump_outlet_region_.empty()) {
      auto pump_outlet_id_list = mesh.getSetEntities(
      pump_outlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  }

  //** ASPIRATIONAL IMPLEMENTATION  -------------------------------
  // // Calculate Pump flowrate
  // // average stage at inlet
  // double avg_pe_inlet = computeAreaWeightedAverage(mesh, pump_inlet_id_list, cv, pe);

  // // average stage at outlet
  // double avg_pe_outlet = max_elev_pumpline_;
  // if (!pump_outlet_region_.empty()) {
  //     avg_pe_outlet = computeAreaWeightedAverage(mesh, pump_outlet_id_list, cv, pe);
  // }

  // double head_diff = std::max(max_elev_pumpline_, avg_pe_outlet) - avg_pe_inlet;
  // double Q = (*Q_pump_)(std::vector<double>{head_diff}); // m^3/s
  //** -----------------------------------------------------------------

  
    // Calculate Pump flowrate
  // average stage at inlet
  double avg_pe_inlet_terms_l[2] = { 0, 0 }; 
  for (auto c : pump_inlet_id_list) {
    avg_pe_inlet_terms_l[0] += cv[0][c] * pe[0][c];
    avg_pe_inlet_terms_l[1] += cv[0][c];
  }
  double avg_pe_inlet_terms_g[2] = { 0, 0 };
  mesh.getComm()->SumAll(avg_pe_inlet_terms_l, avg_pe_inlet_terms_g, 2); 
  double avg_pe_inlet = avg_pe_inlet_terms_g[0] / avg_pe_inlet_terms_g[1]; // area weighted average

  // average stage at outlet
  double avg_pe_outlet = max_elev_pumpline_;
  if (!pump_outlet_region_.empty()) {
      double avg_pe_outlet_terms_l[2] = { 0, 0 }; 
      for (auto c : pump_outlet_id_list) {
        avg_pe_outlet_terms_l[0] += cv[0][c] * pe[0][c];
        avg_pe_outlet_terms_l[1] += cv[0][c];
      }
      double avg_pe_outlet_terms_g[2] = { 0, 0 };
      mesh.getComm()->SumAll(avg_pe_outlet_terms_l, avg_pe_outlet_terms_g, 2); 
      double avg_pe_outlet = avg_pe_outlet_terms_g[0] / avg_pe_outlet_terms_g[1]; 
  }

  double avg_pe_on_off = 0; 
  if (on_off_region_.empty()){
      avg_pe_on_off = avg_pe_inlet;
  } else {
      auto on_off_id_list = mesh.getSetEntities(
      on_off_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      
      double avg_pe_on_off_terms_l[2] = { 0, 0 }; 
      for (auto c : pump_inlet_id_list) {
        avg_pe_on_off_terms_l[0] += cv[0][c] * pe[0][c];
        avg_pe_on_off_terms_l[1] += cv[0][c];
      }
      double avg_pe_on_off_terms_g[2] = { 0, 0 };
      mesh.getComm()->SumAll(avg_pe_on_off_terms_l, avg_pe_on_off_terms_g, 2); 
      avg_pe_on_off = avg_pe_on_off_terms_g[0] / avg_pe_on_off_terms_g[1]; 
  }

  // Determine if the pump should be on or off
  // pump_on = false;
  if (avg_pe_on_off > stage_on_) {
      pump_on = true;
  } else if (avg_pe_on_off < stage_off_) {
      pump_on = false;
  }


  if (pump_on) {
    double head_diff = std::max(max_elev_pumpline_, avg_pe_outlet) - avg_pe_inlet;
    double Q = (*Q_pump_)(std::vector<double>{head_diff}); // m^3/s

    // average ponded depth for sink distribution based on the available water
    double avg_pd_inlet_terms_l[2] = { 0, 0 }; 
    for (auto c : pump_inlet_id_list) {
      avg_pd_inlet_terms_l[0] += cv[0][c] * pd[0][c];
      avg_pd_inlet_terms_l[1] += cv[0][c];
    }
    double avg_pd_inlet_terms_g[2] = { 0, 0 };
    mesh.getComm()->SumAll(avg_pd_inlet_terms_l, avg_pd_inlet_terms_g, 2); 


    // Pump flowrate as sink to inlet cells
    for (auto c : pump_inlet_id_list) {
      if (avg_pd_inlet_terms_g[0] != 0) {
        surf_src[0][c] = - Q * liq_den[0][c] * pd[0][c] / avg_pd_inlet_terms_g[0]; // mol/(m^2 * s)
        // surf_src[0][c] = - Q * liq_den[0][c] / avg_pd_terms_g[1]; 
      }
    }

    // Pump flowrate as source to outlet cells
    if (!pump_outlet_region_.empty()) {
      double sum_cv_l = 0; 
      for (auto c : pump_outlet_id_list) {
        sum_cv_l += cv[0][c];
      }
      double sum_cv_g = 0;
      mesh.getComm()->SumAll(&sum_cv_l, &sum_cv_g, 1);

      for (auto c : pump_outlet_id_list) {
        surf_src[0][c] = Q * liq_den[0][c] / (sum_cv_g); // mol/(m^2 * s)
      }
    }
  }

}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi 