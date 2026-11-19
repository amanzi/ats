/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon  (coonet@ornl.gov)
*/

#include "flux_divergence_evaluator.hh"

namespace Amanzi {
namespace Relations {

FluxDivergenceEvaluator::FluxDivergenceEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // determine the domain
  Key akey = my_keys_.front().first;
  auto tag = my_keys_.front().second;
  Key domain = Keys::getDomain(akey);

  flux_key_ = Keys::readKey(plist, domain, "flux", "water_flux");
  dependencies_.insert(KeyTag{ flux_key_, tag });
}

Teuchos::RCP<Evaluator>
FluxDivergenceEvaluator::Clone() const
{
  return Teuchos::rcp(new FluxDivergenceEvaluator(*this));
}

void
FluxDivergenceEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first, my_keys_.front().second);
  if (my_fac.Mesh() != Teuchos::null) {
    for (auto dep : dependencies_) {
      auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
      fac.SetMesh(my_fac.Mesh())->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    }
  }
}


void
FluxDivergenceEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  const auto& flux_f = *S.Get<CompositeVector>(flux_key_, tag).ViewComponent("face", true);
  const AmanziMesh::Mesh& m = *result[0]->Mesh();
  auto& div_flux = *result[0]->ViewComponent("cell", false);
  div_flux.PutScalar(0.);

  for (AmanziMesh::Entity_ID c = 0; c != div_flux.MyLength(); ++c) {
    auto [cfaces, cfdirs] = m.getCellFacesAndDirections(c);
    for (int i = 0; i != cfaces.size(); ++i) {
      div_flux[0][c] += flux_f[0][cfaces[i]] * cfdirs[i];
    }
  }
}


} // namespace Relations
} // namespace Amanzi
