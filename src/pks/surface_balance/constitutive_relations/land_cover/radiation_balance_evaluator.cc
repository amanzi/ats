/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates a net radiation balance for surface, snow, and canopy.
#include "radiation_balance_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

const std::string RadiationBalanceEvaluator::eval_type = "radiation balance, surface, snow, and canopy";

RadiationBalanceEvaluator::RadiationBalanceEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorModelCVByMaterial<RadiationBalanceModel>(plist)
{}


Teuchos::RCP<Evaluator> RadiationBalanceEvaluator::Clone() const
{
  return Teuchos::rcp(new RadiationBalanceEvaluator(*this));
}

void
RadiationBalanceEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  const auto& deps = models_[0].second->getDependencies();
  Key albedo_key = deps[0].first;
  Key emiss_key = deps[1].first;
  Key area_frac_key = deps[7].first;

  for (const auto& dep : dependencies_) {
    // dependencies on same mesh, but some have two
    if (dep.first == albedo_key || dep.first == emiss_key ||
        dep.first == area_frac_key) {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
        .SetMesh(S.GetMesh(Keys::getDomain(dep.first)))
        ->SetGhosted(false)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 2);
    } else {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
        .SetMesh(S.GetMesh(Keys::getDomain(dep.first)))
        ->SetGhosted(false)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
