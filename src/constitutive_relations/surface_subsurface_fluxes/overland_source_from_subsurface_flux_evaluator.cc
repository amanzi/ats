/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

*/

#include "overland_source_from_subsurface_flux_evaluator.hh"

namespace Amanzi {
namespace Relations {

OverlandSourceFromSubsurfaceFluxEvaluator::OverlandSourceFromSubsurfaceFluxEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_surf_ = Keys::getDomain(my_keys_.front().first);
  domain_sub_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "subsurface");
  Tag tag = my_keys_.front().second;

  flux_key_ = Keys::readKey(plist_, domain_sub_, "flux", "mass_flux");
  dependencies_.insert(KeyTag{ flux_key_, tag });

  // this can be used by both OverlandFlow PK, which uses a volume basis to
  // conserve mass, or OverlandHeadPK, which uses the standard molar basis.
  // If we're using the volume basis, the subsurface's flux (in mol / s) must
  // be divided by molar density to get m^3 / s.
  volume_basis_ = plist_.get<bool>("volume basis", false);
  if (volume_basis_) {
    dens_key_ = Keys::readKey(plist_, domain_sub_, "molar density key", "molar_density_liquid");
    dependencies_.insert(KeyTag{ dens_key_, tag });
  }
}


Teuchos::RCP<Evaluator>
OverlandSourceFromSubsurfaceFluxEvaluator::Clone() const
{
  return Teuchos::rcp(new OverlandSourceFromSubsurfaceFluxEvaluator(*this));
}


void
OverlandSourceFromSubsurfaceFluxEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  S.Require<CompositeVector, CompositeVectorSpace>(my_key, tag, my_key)
    .SetMesh(S.GetMesh(domain_surf_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S.Require<CompositeVector, CompositeVectorSpace>(flux_key_, tag)
    .SetMesh(S.GetMesh(domain_sub_))
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  if (!dens_key_.empty()) {
    S.Require<CompositeVector, CompositeVectorSpace>(dens_key_, tag)
      .SetMesh(S.GetMesh(domain_sub_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }
}


void
OverlandSourceFromSubsurfaceFluxEvaluator::IdentifyFaceAndDirection_(const State& S)
{
  // grab the meshes
  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S.GetMesh(domain_sub_);
  Teuchos::RCP<const AmanziMesh::Mesh> surface = S.GetMesh(domain_surf_);

  // allocate space for face IDs and directions
  int ncells = surface->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  face_and_dirs_ = Teuchos::rcp(new std::vector<FaceDir>(ncells));

  for (int c = 0; c != ncells; ++c) {
    // Get the face on the subsurface mesh corresponding to the cell
    // of the surface mesh.
    AmanziMesh::Entity_ID domain_face = surface->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

    // Get the direction corresponding to that face wrt its only cell.
    // -- get the cell
    auto cells = subsurface->getFaceCells(domain_face, AmanziMesh::Parallel_kind::OWNED);
    AMANZI_ASSERT(cells.size() == 1);

    // -- Get directions
    const auto& [faces, fdirs] = subsurface->getCellFacesAndDirections(cells[0]);
    int index = std::find(faces.begin(), faces.end(), domain_face) - faces.begin();

    // Put (face,dir) into cached data.
    (*face_and_dirs_)[c] = std::make_pair(domain_face, fdirs[index]);
  }
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
OverlandSourceFromSubsurfaceFluxEvaluator::Evaluate_(const State& S,
                                                     const std::vector<CompositeVector*>& result)
{
  if (face_and_dirs_ == Teuchos::null) { IdentifyFaceAndDirection_(S); }
  auto tag = my_keys_.front().second;

  Teuchos::RCP<const AmanziMesh::Mesh> subsurface = S.GetMesh(domain_sub_);
  const Epetra_MultiVector& flux =
    *S.Get<CompositeVector>(flux_key_, tag).ViewComponent("face", false);
  const Epetra_MultiVector& res_v = *result[0]->ViewComponent("cell", false);

  if (volume_basis_) {
    const Epetra_MultiVector& dens =
      *S.Get<CompositeVector>(dens_key_, tag).ViewComponent("cell", false);

    int ncells = result[0]->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      auto cells = subsurface->getFaceCells(
        (*face_and_dirs_)[c].first, AmanziMesh::Parallel_kind::OWNED);
      AMANZI_ASSERT(cells.size() == 1);

      res_v[0][c] =
        flux[0][(*face_and_dirs_)[c].first] * (*face_and_dirs_)[c].second / dens[0][cells[0]];
    }
  } else {
    int ncells = result[0]->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      auto cells = subsurface->getFaceCells(
        (*face_and_dirs_)[c].first, AmanziMesh::Parallel_kind::OWNED);
      AMANZI_ASSERT(cells.size() == 1);

      res_v[0][c] = flux[0][(*face_and_dirs_)[c].first] * (*face_and_dirs_)[c].second;
    }
  }
}

void
OverlandSourceFromSubsurfaceFluxEvaluator::EvaluatePartialDerivative_(
  const State& S,
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
