/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/

/*
  An evaluator for converting the darcy flux to volumetric flux

*/

#include "volumetric_darcy_flux_evaluator.hh"

namespace Amanzi {
namespace Relations {

Volumetric_FluxEvaluator::Volumetric_FluxEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key my_key = my_keys_.front().first;
  auto domain_name = Keys::getDomain(my_key);
  Tag tag = my_keys_.front().second;

  flux_key_ = Keys::readKey(plist_, domain_name, "flux key", "darcy_flux");
  dependencies_.insert(KeyTag{ flux_key_, tag });

  dens_key_ = Keys::readKey(plist_, domain_name, "molar density key", "molar_density_liquid");
  dependencies_.insert(KeyTag{ dens_key_, tag });
}

Teuchos::RCP<Evaluator>
Volumetric_FluxEvaluator::Clone() const
{
  return Teuchos::rcp(new Volumetric_FluxEvaluator(*this));
}


void
Volumetric_FluxEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  auto domain_name = Keys::getDomain(my_key);
  Tag tag = my_keys_.front().second;
  S.Require<CompositeVector, CompositeVectorSpace>(flux_key_, tag)
    .SetMesh(S.GetMesh(domain_name))
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  S.Require<CompositeVector, CompositeVectorSpace>(dens_key_, tag)
    .SetMesh(S.GetMesh(domain_name))
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}

void
Volumetric_FluxEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& darcy_flux =
    *S.Get<CompositeVector>(flux_key_, tag).ViewComponent("face", false);
  const Epetra_MultiVector& molar_density =
    *S.Get<CompositeVector>(dens_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& res_v = *result[0]->ViewComponent("face", false);

  const auto& mesh = *result[0]->Mesh();
  int nfaces_owned = mesh.getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    auto cells = mesh.getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    double n_liq = 0.;
    for (int c = 0; c < cells.size(); c++) n_liq += molar_density[0][c];
    n_liq /= cells.size();
    if (n_liq > 0)
      res_v[0][f] = darcy_flux[0][f] / n_liq;
    else
      res_v[0][f] = 0.;
  }
}


void
Volumetric_FluxEvaluator::EvaluatePartialDerivative_(const State& S,
                                                     const Key& wrt_key,
                                                     const Tag& wrt_tag,
                                                     const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
  // this would require differentiating flux wrt pressure, which we
  // don't do for now.
}

} // namespace Relations
} // namespace Amanzi
