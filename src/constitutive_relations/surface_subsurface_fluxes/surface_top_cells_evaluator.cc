/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

*/

#include "surface_top_cells_evaluator.hh"

namespace Amanzi {
namespace Relations {


SurfaceTopCellsEvaluator::SurfaceTopCellsEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  auto domain_name = Keys::getDomain(my_keys_.front().first);
  auto subsurf_domain = Keys::readDomainHint(plist_, domain_name, "surface", "subsurface");
  auto tag = my_keys_.front().second;

  dependency_key_ =
    Keys::readKey(plist_, subsurf_domain, "subsurface", Keys::getVarName(my_keys_.front().first));
  dependencies_.insert(KeyTag{ dependency_key_, tag });
}


Teuchos::RCP<Evaluator>
SurfaceTopCellsEvaluator::Clone() const
{
  return Teuchos::rcp(new SurfaceTopCellsEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
SurfaceTopCellsEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> sub_vector = S.GetPtr<CompositeVector>(dependency_key_, tag);
  const Epetra_MultiVector& sub_vector_cells = *sub_vector->ViewComponent("cell", false);
  Epetra_MultiVector& result_cells = *result[0]->ViewComponent("cell", false);


  int ncells_surf =
    result[0]->Mesh()->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (unsigned int c = 0; c != ncells_surf; ++c) {
    // get the face on the subsurface mesh
    AmanziMesh::Entity_ID f = result[0]->Mesh()->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

    // get the cell interior to the face
    auto cells = sub_vector->Mesh()->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    AMANZI_ASSERT(cells.size() == 1);

    result_cells[0][c] = sub_vector_cells[0][cells[0]];
  }
}


void
SurfaceTopCellsEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  auto domain_name = Keys::getDomain(my_keys_.front().first);

  CompositeVectorSpace fac;
  fac.SetMesh(S.GetMesh(domain_name)->getParentMesh())->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S, fac);
}


} // namespace Relations
} // namespace Amanzi
