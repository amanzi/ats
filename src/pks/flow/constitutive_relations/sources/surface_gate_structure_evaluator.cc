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
#include "surface_gate_structure_evaluator.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


SurfGateEvaluator::SurfGateEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  pd_key_ = Keys::readKey(plist, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });
  wc_key_ = Keys::readKey(plist, domain_, "water content", "water_content");
  dependencies_.insert(KeyTag{ wc_key_, tag });
  liq_den_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  gate_intake_region_ = plist.get<std::string>("gate intake region");
  storage_region_ = plist.get<std::string>("storage area region");
  // FunctionFactory fac;
  // Q_gate_ = Teuchos::rcp(fac.Create(plist.sublist("function")));
  // gate_func_key_ =
  //   Keys::readKey(*plist_, domain_, "gate function", "gate_function");

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
  

  auto& surf_src= *result[0]->ViewComponent("cell"); // not being reference

  double total = 0.0;
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();

  auto gate_intake_id_list = mesh.getSetEntities(
    gate_intake_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  auto storage_id_list = mesh.getSetEntities(
    storage_region_, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // Calculate the flow rate from gate (from reach into detention pond) 
  double avg_pd_terms_l[2] = { 0, 0 }; 
  for (auto c : gate_intake_id_list) {
    avg_pd_terms_l[0] += cv[0][c] * pd[0][c];
    avg_pd_terms_l[1] += cv[0][c];
  }
  double avg_pd_terms_g[2] = { 0, 0 };
  mesh.getComm()->SumAll(avg_pd_terms_l, avg_pd_terms_g, 2); // move 2 to first place if compiler complains
  double avg_pd = avg_pd_terms_g[0] / avg_pd_terms_g[1];
  double Q = (*Q_gate_)(std::vector<double>{avg_pd}); // m^3/s
  
  // Sink to the reach cells
  double water_avail = 0;
  for (auto c : gate_intake_id_list) {
    water_avail += wc[0][c]; //mols
  }

  double water_demand = Q * dt * liq_den[0][0] ; //mols

  for (auto c : gate_intake_id_list) {
    if (avg_pd_terms_g[0] != 0) {
      surf_src[0][c] = - Q * liq_den[0][c] * pd[0][c] / avg_pd_terms_g[0]; // mol/(m^2 * s)
      // surf_src[0][c] = - Q * liq_den[0][c] / avg_pd_terms_g[1]; 
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