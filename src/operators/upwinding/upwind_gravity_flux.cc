/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

// -----------------------------------------------------------------------------
// ATS
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

UpwindGravityFlux::UpwindGravityFlux(const std::string& pkname,
                                     const Tag& tag,
                                     const Teuchos::RCP<std::vector<WhetStone::Tensor>> K)
  : pkname_(pkname), tag_(tag), K_(K){};

void
UpwindGravityFlux::Update(const CompositeVector& cells,
                          CompositeVector& faces,
                          const State& S,
                          const Teuchos::Ptr<Debugger>& db) const
{
  const auto& g_vec = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  CalculateCoefficientsOnFaces(cells, "cell", g_vec, faces, "face");
};

void
UpwindGravityFlux::Update(const CompositeVector& cells,
                          const std::string cell_component,
                          CompositeVector& faces,
                          const std::string face_component,
                          const State& S,
                          const Teuchos::Ptr<Debugger>& db) const
{
  const auto& g_vec = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  CalculateCoefficientsOnFaces(cells, cell_component, g_vec, faces, face_component);
};


void
UpwindGravityFlux::CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                                const std::string cell_component,
                                                const AmanziGeometry::Point& gravity,
                                                CompositeVector& face_coef,
                                                const std::string face_component) const
{
  double flow_eps = 1.e-10;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.Mesh();

  // initialize the face coefficients
  face_coef.ViewComponent(face_component, true)->PutScalar(0.0);
  if (face_coef.HasComponent("cell")) { face_coef.ViewComponent("cell", true)->PutScalar(1.0); }

  // Note that by scattering, and then looping over all Parallel_kind::ALL cells, we
  // end up getting the correct upwind values in all faces (owned or
  // not) bordering an owned cell.  This is the necessary data for
  // making the local matrices in MFD, so there is no need to
  // communicate the resulting face coeficients.

  // communicate ghosted cells
  cell_coef.ScatterMasterToGhosted(cell_component);

  Epetra_MultiVector& face_coef_v = *face_coef.ViewComponent(face_component, true);
  const Epetra_MultiVector& cell_coef_v = *cell_coef.ViewComponent(cell_component, true);


  for (unsigned int c = 0; c != cell_coef.size(cell_component, true); ++c) {
    const auto& [faces, dirs] = mesh->getCellFacesAndDirections(c);
    AmanziGeometry::Point Kgravity = (*K_)[c] * gravity;

    for (unsigned int n = 0; n != faces.size(); ++n) {
      int f = faces[n];

      const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
      if ((normal * Kgravity) * dirs[n] >= flow_eps) {
        face_coef_v[0][f] = cell_coef_v[0][c];
      } else if (std::abs((normal * Kgravity) * dirs[n]) < flow_eps) {
        face_coef_v[0][f] += cell_coef_v[0][c] / 2.;
      }
    }
  }
};


} // namespace Operators
} // namespace Amanzi
