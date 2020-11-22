/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

ElevationEvaluator::ElevationEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist),
    updated_once_(false), 
    dynamic_mesh_(false) {

  Key domain = Keys::getDomain(plist_.get<std::string>("evaluator name"));
  my_keys_.push_back(Keys::readKey(plist_, domain, "elevation", "elevation"));
  my_keys_.push_back(Keys::readKey(plist_, domain, "slope magnitude", "slope_magnitude"));
  my_keys_.push_back(Keys::readKey(plist_, domain, "aspect", "aspect"));

  // If the mesh changes dynamically (e.g. due to the presence of a deformation
  // pk, then we must make sure that elevation is recomputed every time the 
  // mesh has been deformed. The indicator for the mesh deformation event is the 
  // the deformation field.
  dynamic_mesh_ = plist_.get<bool>("dynamic mesh",false);
  deformation_key_ = Keys::getKey(domain, "deformation");
  if (dynamic_mesh_) dependencies_.insert(deformation_key_);
}

void ElevationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  EvaluateElevationAndSlope_(S, results);

  // If boundary faces are requested, grab the slopes on the internal cell
  Teuchos::Ptr<CompositeVector> slope = results[1];

  if (slope->HasComponent("boundary_face")) {
    const Epetra_Map& vandelay_map = slope->Mesh()->exterior_face_map(false);
    const Epetra_Map& face_map = slope->Mesh()->face_map(false);
    Epetra_MultiVector& slope_bf = *slope->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& slope_c = *slope->ViewComponent("cell",false);

    // calculate boundary face values
    AmanziMesh::Entity_ID_List cells;
    int nbfaces = slope_bf.MyLength();
    for (int bf=0; bf!=nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      slope->Mesh()->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      slope_bf[0][bf] = slope_c[0][cells[0]];
    }
  }
}

// This is hopefully never called?
void ElevationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  AMANZI_ASSERT(0);
}

// Custom EnsureCompatibility forces this to be updated once.
bool ElevationEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
                                         Key request) {
  bool changed = SecondaryVariablesFieldEvaluator::HasFieldChanged(S,request);
  if (!updated_once_) {
    UpdateField_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}


void ElevationEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<CompositeVectorSpace> master_fac;
  for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
       my_key!=my_keys_.end(); ++my_key) {
    // Ensure my field exists, and claim ownership.
    Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(*my_key, *my_key);

    // Check plist for vis or checkpointing control.
    bool io_my_key = plist_.get<bool>(std::string("visualize ")+*my_key, true);
    S->GetField(*my_key, *my_key)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+*my_key, false);
    S->GetField(*my_key, *my_key)->set_io_checkpoint(checkpoint_my_key);

    // Select a "master factory" to ensure commonality of all of my factories.
    if (my_fac->Mesh() != Teuchos::null && master_fac == Teuchos::null) {
      master_fac = my_fac;
    }
  }

  if (master_fac == Teuchos::null) {
    // No requirements have been set, so we'll have to hope they get set by an
    // evaluator that depends upon this evaluator.
  } else {
    for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
         my_key!=my_keys_.end(); ++my_key) {
      Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(*my_key, *my_key);

      // Cannot just Update() the factory, since the locations may differ, but
      // at least we can ensure the mesh exists.
      my_fac->SetMesh(master_fac->Mesh());
      my_fac->AddComponent("cell", AmanziMesh::CELL, 1);
    }
  }
};



} //namespace
} //namespace
