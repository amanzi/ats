/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

/*
This evaluator models gate structure inside 2D flow area. The gate curve is kept general and can be provided by the user.
Gate structure can be used to move water between two canals, two storage areas or canal to storage area.
*/

#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "surface_gate_structure_evaluator.hh"
#include "FunctionFactory.hh"


namespace Amanzi {
namespace Flow {
namespace Relations {

// Helper function to compute area-weighted average
inline double
computeAreaWeightedAverage(const AmanziMesh::Mesh& mesh,
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


SurfGateEvaluator::SurfGateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
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

  gate_intake_region_ = plist.get<std::string>("gate intake region");
  storage_region_ = plist.get<std::string>("storage area region");
  is_ponded_depth_function_ = plist.get<bool>("is ponded depth function", false);
  stage_close_ = plist.get<double>("gate close stage", 1e6);

  Teuchos::ParameterList& gate_func = plist.sublist("function");
  FunctionFactory fac;
  Q_gate_ = Teuchos::rcp(fac.Create(gate_func));
}


// Required methods from SecondaryVariableFieldEvaluator
void
SurfGateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
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

  auto gate_intake_id_list = mesh.getSetEntities(
    gate_intake_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  auto storage_id_list = mesh.getSetEntities(
    storage_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // Calculate the flow rate from gate (from reach into detention pond)
  double Q;
  double avg_pe_storage = computeAreaWeightedAverage(mesh, storage_id_list, cv, pe);

  if (is_ponded_depth_function_) {
    double avg_pd_inlet = computeAreaWeightedAverage(mesh, gate_intake_id_list, cv, pd);
    // test that for zero ponded depth, the function returns zero flow
    double Q0 = (*Q_gate_)(std::vector<double>{0});
    AMANZI_ASSERT(Q0 == 0);
    Q = (*Q_gate_)(std::vector<double>{avg_pd_inlet}); // m^3/s
  } else {
    double avg_pe_inlet = computeAreaWeightedAverage(mesh, gate_intake_id_list, cv, pe);
    double avg_elev_inlet = computeAreaWeightedAverage(mesh, gate_intake_id_list, cv, elev);
    double Q0 = (*Q_gate_)(std::vector<double>{avg_elev_inlet});
    AMANZI_ASSERT(Q0 == 0);
    Q = (*Q_gate_)(std::vector<double>{avg_pe_inlet});
  }

  // to overcome "bang-bang" behavior, we use a smooth transition
  double delta_ = 0.01; // should we make it available to the user?
  double alpha = 0.5 * (1.0 - std::tanh((avg_pe_storage - stage_close_) / delta_));
  Q = alpha * Q;

  // Sink to the reach cells
  double sum_wc_l = 0;
  for (auto c : gate_intake_id_list) {
    sum_wc_l += wc[0][c];
  }
  double sum_wc_g = 0;
  mesh.getComm()->SumAll(&sum_wc_l, &sum_wc_g, 1);

  for (auto c : gate_intake_id_list) {
    if (sum_wc_g > 0) {
      surf_src[0][c] = -Q * liq_den[0][c] * wc[0][c] / (sum_wc_g * cv[0][c]); // mol/(m^2 * s)
    }
  }

  // Source to the detention pond cells
  double sum_cv_l = 0;
  for (auto c : storage_id_list) {
    sum_cv_l += cv[0][c];
  }
  double sum_cv_g = 0;
  mesh.getComm()->SumAll(&sum_cv_l, &sum_cv_g, 1);

  for (auto c : storage_id_list) {
    surf_src[0][c] = Q * liq_den[0][c] / (sum_cv_g); // mol/(m^2 * s)
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi