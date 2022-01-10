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
    SecondaryVariableFieldEvaluator(plist) {
  auto domain = Keys::getDomain(my_key_);
  std::string surf_domain;
  if (Keys::getDomain(my_key_).empty()) {
    surf_domain = "surface";
  } else {
    surf_domain = Key("surface_")+domain;
  }
  surf_domain = plist_.get<std::string>("surface domain name", surf_domain);
  dependency_key_ = Keys::readKey(plist_, surf_domain, "surface", Keys::getVarName(my_key_));
  dependencies_.insert(dependency_key_);

  negate_ = plist_.get<bool>("negate", false);
}

TopCellsSurfaceEvaluator::TopCellsSurfaceEvaluator(const TopCellsSurfaceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    negate_(other.negate_),
    dependency_key_(other.dependency_key_) {}

Teuchos::RCP<FieldEvaluator>
TopCellsSurfaceEvaluator::Clone() const {
  return Teuchos::rcp(new TopCellsSurfaceEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
TopCellsSurfaceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> surf_vector = S->GetPtr<CompositeVector>(dependency_key_);
  const Epetra_MultiVector& surf_vector_cells = *surf_vector->ViewComponent("cell",false);
  Epetra_MultiVector& result_cells = *result->ViewComponent("cell",false);


  int ncells_surf = surf_vector->Mesh()->num_entities(AmanziMesh::CELL,
          AmanziMesh::Parallel_type::OWNED);
  for (unsigned int c=0; c!=ncells_surf; ++c) {
    // get the face on the subsurface mesh
    AmanziMesh::Entity_ID f = surf_vector->Mesh()->entity_get_parent(AmanziMesh::CELL, c);

    // get the cell interior to the face
    AmanziMesh::Entity_ID_List cells;
    result->Mesh()->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);

    result_cells[0][cells[0]] = surf_vector_cells[0][c];
  }
  if (negate_) result->Scale(-1);
}


void
TopCellsSurfaceEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Ensure my field exists.  Requirements should be already set.  Claim ownership.
  AMANZI_ASSERT(my_key_ != std::string(""));
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT,  my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
}



} //namespace
} //namespace

