/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#include "Tensor.hh"
#include "CompositeVector.hh"
#include "State.hh"
#include "upwind_gravity_flux.hh"

namespace Amanzi {
namespace Operators {

UpwindGravityFlux:: UpwindGravityFlux(const std::string& pkname,
        const Tag& tag,
        const Teuchos::RCP<std::vector<WhetStone::Tensor> > K)
  : pkname_(pkname),
    tag_(tag),
    K_(K) {};


void
UpwindGravityFlux::Update(const CompositeVector& cells,
                          CompositeVector& faces,
                          const State& S,
                          const Teuchos::Ptr<Debugger>& db) const
{
  const auto& g_vec = S.Get<AmanziGeometry::Point>("gravity");
  CalculateCoefficientsOnFaces(cells, g_vec, faces);
};


void UpwindGravityFlux::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const AmanziGeometry::Point& gravity,
        CompositeVector& face_coef) const
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  double flow_eps = 1.e-10;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.Mesh();

  // initialize the face coefficients
  face_coef.ViewComponent("face",true)->PutScalar(0.0);
  if (face_coef.HasComponent("cell")) {
    face_coef.ViewComponent("cell",true)->PutScalar(1.0);
  }

  // Note that by scattering, and then looping over all Parallel_type::ALL cells, we
  // end up getting the correct upwind values in all faces (owned or
  // not) bordering an owned cell.  This is the necessary data for
  // making the local matrices in MFD, so there is no need to
  // communicate the resulting face coeficients.

  // communicate ghosted cells
  cell_coef.ScatterMasterToGhosted("cell");

  Epetra_MultiVector& face_coef_v = *face_coef.ViewComponent("face",true);
  const Epetra_MultiVector& cell_coef_v = *cell_coef.ViewComponent("cell",true);


  for (unsigned int c=0; c!=cell_coef.size("cell", true); ++c) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    AmanziGeometry::Point Kgravity = (*K_)[c] * gravity;

    for (unsigned int n=0; n!=faces.size(); ++n) {
      int f = faces[n];

      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      if ((normal * Kgravity) * dirs[n] >= flow_eps) {
        face_coef_v[0][f] = cell_coef_v[0][c];
      } else if (std::abs((normal * Kgravity) * dirs[n]) < flow_eps) {
        face_coef_v[0][f] += cell_coef_v[0][c] / 2.;
      }
    }
  }
};


} //namespace
} //namespace
