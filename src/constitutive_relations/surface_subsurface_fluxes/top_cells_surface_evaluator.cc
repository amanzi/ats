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

#include "top_cells_surface_evaluator.hh"

namespace Amanzi {
namespace Relations {


TopCellsSurfaceEvaluator::TopCellsSurfaceEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key my_key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_key);
  domain_surf_ = Keys::readDomainHint(plist_, domain, "subsurface", "surface");
  dependency_key_ = Keys::readKey(plist_, domain_surf_, "surface", Keys::getVarName(my_key));
  dependencies_.insert(KeyTag{ dependency_key_, tag });

  negate_ = plist_.get<bool>("negate", false);
}

Teuchos::RCP<Evaluator>
TopCellsSurfaceEvaluator::Clone() const
{
  return Teuchos::rcp(new TopCellsSurfaceEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotype
void
TopCellsSurfaceEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  auto surf_vector = S.GetPtr<CompositeVector>(dependency_key_, tag);
  const Epetra_MultiVector& surf_vector_cells = *surf_vector->ViewComponent("cell", false);
  Epetra_MultiVector& result_cells = *result[0]->ViewComponent("cell", false);

  int ncells_surf =
    surf_vector->Mesh()->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (unsigned int c = 0; c != ncells_surf; ++c) {
    // get the face on the subsurface mesh
    AmanziMesh::Entity_ID f = surf_vector->Mesh()->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

    // get the cell interior to the face
    auto cells = result[0]->Mesh()->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    AMANZI_ASSERT(cells.size() == 1);

    result_cells[0][cells[0]] = surf_vector_cells[0][c];
  }
  if (negate_) result[0]->Scale(-1);
}


void
TopCellsSurfaceEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  CompositeVectorSpace fac;
  fac.SetMesh(S.GetMesh(domain_surf_))->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S, fac);
}


} // namespace Relations
} // namespace Amanzi
