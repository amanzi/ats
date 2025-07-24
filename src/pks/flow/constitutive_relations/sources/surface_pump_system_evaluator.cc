/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)
           Ethan Coon (coonet@ornl.gov)

*/

/*
This evaluator models stage-based pump station model inside 2D flow area.
Pump stations can be used to move water between any combination of river reaches, storage areas or catchment regions.
Based on pump on-off conditions and pump-operation curve, water is moved from pump-inlet to -outlet region instantly.
*/

#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "surface_pump_system_evaluator.hh"
#include "FunctionFactory.hh"
#include <Epetra_MultiVector.h>


namespace Amanzi {
namespace Flow {
namespace Relations {

static const double NaN = std::numeric_limits<double>::signaling_NaN();

// Helper function to compute area-weighted average
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


SurfPumpEvaluator::SurfPumpEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  auto domain = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  // dependencies
  cv_key_ = Keys::readKey(plist, domain, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  pd_key_ = Keys::readKey(plist, domain, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });

  pe_key_ = Keys::readKey(plist, domain, "potential", "pres_elev");
  dependencies_.insert(KeyTag{ pe_key_, tag });

  elev_key_ = Keys::readKey(plist, domain, "elevation", "elevation");
  dependencies_.insert(KeyTag{ elev_key_, tag });

  wc_key_ = Keys::readKey(plist, domain, "water content", "water_content");
  dependencies_.insert(KeyTag{ wc_key_, tag });

  liq_den_key_ = Keys::readKey(plist, domain, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ liq_den_key_, tag });

  // need an extra flag stored in state to indiciate, from
  // timestep-to-timestep, that the pump is on or off.
  pump_on_key_ = Keys::readKey(
    plist, domain, "pump on", Keys::getVarName(my_keys_.front().first) + "_pump_on_flag");

  pump_outlet_region_ = plist.get<std::string>("pump outlet region", "");
  pump_inlet_region_ = plist.get<std::string>("pump inlet region");
  on_off_region_ = plist.get<std::string>("on off reference region", "");

  max_elev_pumpline_ = plist.get<double>("maximum pumpline elevation", NaN);
  stage_on_ = plist.get<double>("pump start at stage", NaN);
  stage_off_ = plist.get<double>("pump stop at stage", NaN);

  Teuchos::ParameterList& pump_func = plist.sublist("function");
  FunctionFactory fac;
  Q_pump_ = Teuchos::rcp(fac.Create(pump_func));
}


void
SurfPumpEvaluator::EnsureCompatibility_Structure_(State& S)
{
  auto tag = my_keys_.front().second;
  S.Require<int>(pump_on_key_, tag, my_keys_.front().first);

  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Structure_(S);
}


void
SurfPumpEvaluator::Update_(State& S)
{
  // vector of pointers to results
  std::vector<CompositeVector*> results;
  for (const auto& keytag : my_keys_) {
    results.push_back(&S.GetW<CompositeVector>(keytag.first, keytag.second, keytag.first));
  }
  // call the evaluate method
  int& pump_on = S.GetW<int>(pump_on_key_, my_keys_.front().second, my_keys_.front().first);
  Evaluate_(S, results, pump_on);
}


// Required methods from SecondaryVariableFieldEvaluator
void
SurfPumpEvaluator::Evaluate_(const State& S,
                             const std::vector<CompositeVector*>& result,
                             int& pump_on)
{
  Tag tag = my_keys_.front().second;

  double dt = S.Get<double>("dt", tag);
  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& pd = *S.Get<CompositeVector>(pd_key_, tag).ViewComponent("cell", false);
  const auto& liq_den = *S.Get<CompositeVector>(liq_den_key_, tag).ViewComponent("cell", false);
  const auto& wc = *S.Get<CompositeVector>(wc_key_, tag).ViewComponent("cell", false);
  const auto& pe = *S.Get<CompositeVector>(pe_key_, tag).ViewComponent("cell", false);
  const auto& elev = *S.Get<CompositeVector>(elev_key_, tag).ViewComponent("cell", false);

  auto& surf_src = *result[0]->ViewComponent("cell");

  double total = 0.0;
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  surf_src.PutScalar(0.); // initializing with zero

  // collect gids in relevant regions
  auto pump_inlet_id_list = mesh.getSetEntities(
    pump_inlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  Amanzi::AmanziMesh::MeshCache<Amanzi::MemSpace_kind::HOST>::cEntity_ID_View pump_outlet_id_list;
  if (!pump_outlet_region_.empty()) {
    auto pump_outlet_id_list = mesh.getSetEntities(
      pump_outlet_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  }

  // stage_on value should be greater than elevation at inlet or inlet is the reference region
  double avg_elev_inlet;
  if (on_off_region_.empty()) {
    avg_elev_inlet = computeAreaWeightedAverage(mesh, pump_inlet_id_list, cv, elev);
    AMANZI_ASSERT(stage_on_ > avg_elev_inlet);
  }

  // Calculate Stage Difference and Pump Status
  // average stage at inlet
  double avg_pe_inlet = computeAreaWeightedAverage(mesh, pump_inlet_id_list, cv, pe);
  double avg_pd_inlet = computeAreaWeightedAverage(mesh, pump_inlet_id_list, cv, pd);

  // average stage at outlet if outlet region provided
  double avg_pe_outlet;
  if (!pump_outlet_region_.empty()) {
    avg_pe_outlet = computeAreaWeightedAverage(mesh, pump_outlet_id_list, cv, pe);
  } else {
    // if the pump is at the domain boundary pumping water out of the domain,
    // need to provide max pumpline elevation
    avg_pe_outlet = max_elev_pumpline_;
  }

  double avg_pe_on_off = 0;
  if (on_off_region_.empty()) {
    avg_pe_on_off = avg_pe_inlet;
  } else { // if on off region is separate from the inlet region
    auto on_off_id_list = mesh.getSetEntities(
      on_off_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    avg_pe_on_off = computeAreaWeightedAverage(mesh, on_off_id_list, cv, pe);
  }

  // Determine if the pump should be on or off
  if (avg_pe_on_off > stage_on_) {
    pump_on = 1;
  } else if (avg_pe_on_off < stage_off_) {
    pump_on = 0;
  } else {
    // pump is left at the same state as last timestep
    // pass
  }

  // Calculate pumpflow rate and distribute sources and sinks
  double Q_max;
  double Q;
  if (pump_on) {
    double head_diff = std::max(max_elev_pumpline_, avg_pe_outlet) - avg_pe_inlet;
    Q_max = (*Q_pump_)(std::vector<double>{ head_diff }); // m^3/s
    Q = Q_max;

    // sink distribution proportional to available water
    double sum_wc_l = 0;
    for (auto c : pump_inlet_id_list) {
      sum_wc_l += wc[0][c];
    }
    double sum_wc_g = 0;
    mesh.getComm()->SumAll(&sum_wc_l, &sum_wc_g, 1);

    // Pump flowrate as sink to inlet cells
    for (auto c : pump_inlet_id_list) {
      if (sum_wc_g != 0) {
        surf_src[0][c] = -Q * liq_den[0][c] * wc[0][c] / (sum_wc_g * cv[0][c]); // mol/(m^2 * s)
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

  // double-check mass conservation
  double mass_balance = 0.0;
  if (!pump_outlet_region_.empty()) {
    surf_src.Dot(cv, &mass_balance);
    // Check if the mass balance is zero
    AMANZI_ASSERT(std::abs(mass_balance) < 1.e-10);
  } else { // Debug code
    for (auto c : pump_inlet_id_list) {
      mass_balance += surf_src[0][c] * cv[0][c] / liq_den[0][c]; // m^3/s
    }
    double mass_balance_g = 0;
    mesh.getComm()->SumAll(&mass_balance, &mass_balance_g, 1);
    AMANZI_ASSERT(std::abs(mass_balance_g + Q) < 1.e-10);
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
