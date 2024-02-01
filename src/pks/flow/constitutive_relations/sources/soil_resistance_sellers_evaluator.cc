/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*
  Evaluates the soil resistance at top cells through the Sellers method
  and assign them to surface cells.
*/


#include "MeshAlgorithms.hh"
#include "soil_resistance_sellers_evaluator.hh"

namespace Amanzi {
namespace Flow {

SoilResistanceSellersEvaluator::SoilResistanceSellersEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  std::string domain_surf = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  Key domain_ss = Keys::readDomainHint(plist, domain_surf, "surface", "subsurface");

  // my dependencies
  sat_key_ = Keys::readKey(plist_, domain_ss, "liquid saturation", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_key_, tag });
}


Teuchos::RCP<Evaluator>
SoilResistanceSellersEvaluator::Clone() const
{
  return Teuchos::rcp(new SoilResistanceSellersEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
SoilResistanceSellersEvaluator::Evaluate_(const State& S,
                                          const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> sat = S.GetPtr<CompositeVector>(sat_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = sat->Mesh();
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = result[0]->Mesh();

  // evaluate the model
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell"); // partition on cell only
    const Epetra_MultiVector& sat_v = *(sat->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int sc = 0; sc != count; ++sc) {
      int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
      result_v[0][sc] = std::exp(8.206 - 4.255 * sat_v[0][c]);
    }
  }
}


void
SoilResistanceSellersEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> sat = S.GetPtr<CompositeVector>(sat_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = sat->Mesh();
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = result[0]->Mesh();

  if (wrt_key == sat_key_) {
    // evaluate the model
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      AMANZI_ASSERT(*comp == "cell"); // partition on cell only
      const Epetra_MultiVector& sat_v = *(sat->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int sc = 0; sc != count; ++sc) {
        int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        result_v[0][sc] = -4.255 * std::exp(8.206 - 4.255 * sat_v[0][c]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}


void
SoilResistanceSellersEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  const auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                     my_keys_.front().second);
  if (fac.Mesh() != Teuchos::null) {
    CompositeVectorSpace dep_fac;
    dep_fac.SetMesh(fac.Mesh()->getParentMesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    for (const auto& dep : dependencies_) {
      if (Keys::getDomain(dep.first) == Keys::getDomain(my_keys_.front().first)) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(dep_fac);
      }
    }
  }
}


} // namespace Flow
} // namespace Amanzi
