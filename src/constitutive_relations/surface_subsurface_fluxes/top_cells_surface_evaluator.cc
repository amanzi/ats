/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "top_cells_surface_evaluator.hh"

namespace Amanzi {
namespace Relations {


TopCellsSurfaceEvaluator::TopCellsSurfaceEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  Key my_key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_key);
  Key surf_domain = Keys::readDomainHint(plist_, domain, "subsurface", "surface");
  dependency_key_ = Keys::readKey(plist_, surf_domain, "surface", Keys::getVarName(my_key));
  dependencies_.insert(KeyTag{dependency_key_, tag});

  negate_ = plist_.get<bool>("negate", false);
}

Teuchos::RCP<Evaluator>
TopCellsSurfaceEvaluator::Clone() const {
  return Teuchos::rcp(new TopCellsSurfaceEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotype
void
TopCellsSurfaceEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  auto surf_vector = S.GetPtr<CompositeVector>(dependency_key_);
  const Epetra_MultiVector& surf_vector_cells = *surf_vector->ViewComponent("cell",false);
  Epetra_MultiVector& result_cells = *result[0]->ViewComponent("cell",false);

  int ncells_surf = surf_vector->Mesh()->num_entities(AmanziMesh::CELL,
          AmanziMesh::Parallel_type::OWNED);
  for (unsigned int c=0; c!=ncells_surf; ++c) {
    // get the face on the subsurface mesh
    AmanziMesh::Entity_ID f = surf_vector->Mesh()->entity_get_parent(AmanziMesh::CELL, c);

    // get the cell interior to the face
    AmanziMesh::Entity_ID_List cells;
    result[0]->Mesh()->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);

    result_cells[0][cells[0]] = surf_vector_cells[0][c];
  }
  if (negate_) result[0]->Scale(-1);
}


void
TopCellsSurfaceEvaluator::EnsureCompatibility(State& S)
{
  // Ensure my field exists.  Requirements should be already set.  Claim ownership.
  Key my_key = my_keys_.front().first;
  S.Require<CompositeVector,CompositeVectorSpace>(my_key,
          my_keys_.front().second,  my_key)
    .SetMesh(S.GetMesh(Keys::getDomain(my_key)))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  for (const auto& dep : dependencies_) {
    S.Require<CompositeVector,CompositeVectorSpace>(dep.first, dep.second)
      .SetMesh(S.GetMesh(Keys::getDomain(dep.first)))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.RequireEvaluator(dep.first, dep.second).EnsureCompatibility(S);
  }

  // check plist for vis or checkpointing control
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Flags_(S);
}



} //namespace
} //namespace

