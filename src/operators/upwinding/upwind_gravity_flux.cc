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
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  double flow_eps = 1.e-10;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.getMesh();

  // initialize the face coefficients
  face_coef.getComponent(face_component, true)->putScalar(0.0);
  if (face_coef.hasComponent("cell")) {
    face_coef.getComponent("cell", true)->putScalar(1.0);
  }

  // Note that by scattering, and then looping over all Parallel_kind::ALL cells, we
  // end up getting the correct upwind values in all faces (owned or
  // not) bordering an owned cell.  This is the necessary data for
  // making the local matrices in MFD, so there is no need to
  // communicate the resulting face coeficients.

  // communicate ghosted cells
  cell_coef.scatterMasterToGhosted(cell_component);

  {
    auto face_coef_v = face_coef.viewComponent(face_component, true);
    const auto cell_coef_v = cell_coef.viewComponent(cell_component, true);

    int nfaces_local = face_coef_v.extent(0);
    // NOTE -- this needs a new algorithm -- not continuous
    AMANZI_ASSERT(false);
    // Kokkos::parallel_for("upwind_gravity_flux", nfaces_local,
    //                      KOKKOS_LAMBDA(const int& f) {
    //                        auto fcells = mesh->getFaceCells(f);
    //                        int c0 = fcells(0);
    //                        int orientation = 0;
    //                        auto normal = mesh->getFaceNormal(f, c0, &orientation);

    //                        AmanziGeometry::Point Kgravity = (*K_)[c0] * gravity;

    //                        auto grav_flux = normal * Kgravity;
    //                        if (grav_flux >= flow_eps) {
    //                          face_coef_v[0][f] = cell_coef_v[0][c0];
    //                        } else if (grav_flux ) {
    //                          face_coef_v[0][f] += cell_coef_v[0][c] / 2.;
    //   });
  }
};


} // namespace Operators
} // namespace Amanzi
