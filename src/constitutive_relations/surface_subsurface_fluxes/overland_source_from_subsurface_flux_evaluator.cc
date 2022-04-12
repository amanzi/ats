/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_source_from_subsurface_flux_evaluator.hh"

namespace Amanzi {
namespace Relations {

OverlandSourceFromSubsurfaceFluxEvaluator::OverlandSourceFromSubsurfaceFluxEvaluator(
        Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("source key", "overland_source_from_subsurface");
  }

  flux_key_ = plist_.get<std::string>("flux key", "water_flux");

  // since we cannot have flux as a dependency (it has no model), we have to
  // use pressure as a proxy.
  Key pres_key = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pres_key);

  // this can be used by both OverlandFlow PK, which uses a volume basis to
  // conserve mass, or OverlandHeadPK, which uses the standard molar basis.
  // If we're using the volume basis, the subsurface's flux (in mol / s) must
  // be divided by molar density to get m^3 / s.
  volume_basis_ = plist_.get<bool>("volume basis", false);
  if (volume_basis_) {
    dens_key_ = plist_.get<std::string>("molar density key", "molar_density_liquid");
    dependencies_.insert(dens_key_);
  }

  surface_mesh_key_ = plist_.get<std::string>("surface mesh key", "surface");
  subsurface_mesh_key_ = plist_.get<std::string>("subsurface mesh key", "domain");
}

OverlandSourceFromSubsurfaceFluxEvaluator::OverlandSourceFromSubsurfaceFluxEvaluator(
        const OverlandSourceFromSubsurfaceFluxEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    flux_key_(other.flux_key_),
    dens_key_(other.dens_key_),
    surface_mesh_key_(other.surface_mesh_key_),
    subsurface_mesh_key_(other.subsurface_mesh_key_),
    face_and_dirs_(other.face_and_dirs_),
    volume_basis_(other.volume_basis_) {}

Teuchos::RCP<FieldEvaluator> OverlandSourceFromSubsurfaceFluxEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandSourceFromSubsurfaceFluxEvaluator(*this));
}


void OverlandSourceFromSubsurfaceFluxEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // for now just passing... might do something later here?
  S->RequireField(my_key_, my_key_);
  S->RequireField(flux_key_);
}


void OverlandSourceFromSubsurfaceFluxEvaluator::IdentifyFaceAndDirection_(
        const Teuchos::Ptr<State>& S) {
  // grab the meshes
  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S->GetMesh(subsurface_mesh_key_);
  Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh(surface_mesh_key_);

  // allocate space for face IDs and directions
  int ncells = surface->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  face_and_dirs_ = Teuchos::rcp(new std::vector<FaceDir>(ncells));

  for (int c=0; c!=ncells; ++c) {
    // Get the face on the subsurface mesh corresponding to the cell
    // of the surface mesh.
    AmanziMesh::Entity_ID domain_face =
      surface->entity_get_parent(AmanziMesh::CELL, c);

    // Get the direction corresponding to that face wrt its only cell.
    // -- get the cell
    AmanziMesh::Entity_ID_List cells;
    subsurface->face_get_cells(domain_face, AmanziMesh::Parallel_type::OWNED, &cells);
    AMANZI_ASSERT(cells.size() == 1);

    // -- Get directions
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> fdirs;
    subsurface->cell_get_faces_and_dirs(cells[0], &faces, &fdirs);
    int index = std::find(faces.begin(), faces.end(), domain_face) - faces.begin();

    // Put (face,dir) into cached data.
    (*face_and_dirs_)[c] = std::make_pair(domain_face, fdirs[index]);
  }
}

// Required methods from SecondaryVariableFieldEvaluator
void OverlandSourceFromSubsurfaceFluxEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  if (face_and_dirs_ == Teuchos::null) {
    IdentifyFaceAndDirection_(S);
  }

  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S->GetMesh(subsurface_mesh_key_);
  const Epetra_MultiVector& flux = *S->GetFieldData(flux_key_)->ViewComponent("face",false);
  const Epetra_MultiVector& res_v = *result->ViewComponent("cell",false);

  if (volume_basis_) {
    const Epetra_MultiVector& dens = *S->GetFieldData(dens_key_)->ViewComponent("cell",false);

    int ncells = result->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      AmanziMesh::Entity_ID_List cells;
      subsurface->face_get_cells((*face_and_dirs_)[c].first, AmanziMesh::Parallel_type::OWNED, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      res_v[0][c] = flux[0][(*face_and_dirs_)[c].first] * (*face_and_dirs_)[c].second
          / dens[0][cells[0]];
    }
  } else {
    int ncells = result->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      AmanziMesh::Entity_ID_List cells;
      subsurface->face_get_cells((*face_and_dirs_)[c].first, AmanziMesh::Parallel_type::OWNED, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      res_v[0][c] = flux[0][(*face_and_dirs_)[c].first] * (*face_and_dirs_)[c].second;
    }
  }
}

void OverlandSourceFromSubsurfaceFluxEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(0);
  // this would require differentiating flux wrt pressure, which we
  // don't do for now.
}

} // namespace
} // namespace
