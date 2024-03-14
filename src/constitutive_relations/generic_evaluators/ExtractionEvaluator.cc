/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Extracts a field on one mesh from a field on a superset of that mesh using
parent entities.

*/

#include "ExtractionEvaluator.hh"

namespace Amanzi {
namespace Relations {


ExtractionEvaluator::ExtractionEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key my_key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  parent_domain_ = Keys::readDomain(plist_, "parent");
  dependency_key_ = Keys::readKey(plist_, parent_domain_, "extracted", Keys::getVarName(my_key));
  dependencies_.insert(KeyTag{ dependency_key_, tag });
}

Teuchos::RCP<Evaluator>
ExtractionEvaluator::Clone() const
{
  return Teuchos::rcp(new ExtractionEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotype
void
ExtractionEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  auto& parent_vector = S.Get<CompositeVector>(dependency_key_, tag);

  for (const auto& comp : *result[0]) {
    const Epetra_MultiVector& parent_vector_c = *parent_vector.ViewComponent(comp, false);
    Epetra_MultiVector& result_c = *result[0]->ViewComponent(comp, false);

    AmanziMesh::Entity_kind entity = result[0]->Location(comp);
    auto mesh = result[0]->Mesh();

    for (int j = 0; j != result_c.NumVectors(); ++j) {
      for (int i = 0; i != result_c.MyLength(); ++i) {
        result_c[j][i] = parent_vector_c[j][mesh->getEntityParent(entity, i)];
      }
    }
  }
}


void
ExtractionEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag my_tag = my_keys_.front().second;
  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key, my_tag);

  CompositeVectorSpace fac;
  fac.SetMesh(S.GetMesh(parent_domain_));
  for (auto& comp : my_fac) {
    fac.AddComponent(comp, my_fac.Location(comp), my_fac.NumVectors(comp));
  }
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S, fac);
}


} // namespace Relations
} // namespace Amanzi
