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
// faces. Gives priority to the harmonic average when feasible.
// -----------------------------------------------------------------------------

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "State.hh"
#include "Debugger.hh"
#include "VerboseObject.hh"
#include "upwind_flux_harmonic_mean.hh"
#include "Epetra_IntVector.h"

namespace Amanzi {
namespace Operators {

UpwindFluxHarmonicMean::UpwindFluxHarmonicMean(const std::string& pkname,
                                               const Tag& tag,
                                               const Key& flux,
                                               double flux_epsilon)
  : pkname_(pkname), tag_(tag), flux_(flux), flux_eps_(flux_epsilon)
{}


void
UpwindFluxHarmonicMean::Update(const CompositeVector& cell_coef,
                               CompositeVector& face_coef,
                               const State& S,
                               const Teuchos::Ptr<Debugger>& db) const
{
  const CompositeVector& flux = S.Get<CompositeVector>(flux_, tag_);
  CalculateCoefficientsOnFaces(cell_coef, flux, face_coef, db);
};


void
UpwindFluxHarmonicMean::CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                                     const CompositeVector& flux,
                                                     CompositeVector& face_coef,
                                                     const Teuchos::Ptr<Debugger>& db) const
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.Mesh();

  // initialize the face coefficients
  if (face_coef.HasComponent("cell")) { face_coef.ViewComponent("cell", true)->PutScalar(1.0); }

  // communicate needed ghost values
  cell_coef.ScatterMasterToGhosted("cell");

  // pull out vectors
  const Epetra_MultiVector& flux_v = *flux.ViewComponent("face", false);
  Epetra_MultiVector& coef_faces = *face_coef.ViewComponent("face", false);
  const Epetra_MultiVector& coef_cells = *cell_coef.ViewComponent("cell", true);

  // Identify upwind/downwind cells for each local face.  Note upwind/downwind
  // may be a ghost cell.
  Epetra_IntVector upwind_cell(*face_coef.ComponentMap("face", true));
  upwind_cell.PutValue(-1);
  Epetra_IntVector downwind_cell(*face_coef.ComponentMap("face", true));
  downwind_cell.PutValue(-1);

  int nfaces_local = flux.size("face", false);

  int ncells = cell_coef.size("cell", true);
  for (int c = 0; c != ncells; ++c) {
    const auto& [faces, fdirs] = mesh->getCellFacesAndDirections(c);

    for (unsigned int n = 0; n != faces.size(); ++n) {
      int f = faces[n];

      if (f < nfaces_local) {
        if (flux_v[0][f] * fdirs[n] > 0) {
          upwind_cell[f] = c;
        } else if (flux_v[0][f] * fdirs[n] < 0) {
          downwind_cell[f] = c;
        } else {
          // We don't care, but we have to get one into upwind and the other
          // into downwind.
          if (upwind_cell[f] == -1) {
            upwind_cell[f] = c;
          } else {
            downwind_cell[f] = c;
          }
        }
      }
    }
  }

  // Determine the face coefficient of local faces.
  // These parameters may be key to a smooth convergence rate near zero flux.
  double coefs[2];

  int nfaces = face_coef.size("face", false);
  for (int f = 0; f != nfaces; ++f) {
    int uw = upwind_cell[f];
    int dw = downwind_cell[f];
    AMANZI_ASSERT(!((uw == -1) && (dw == -1)));

    // uw coef
    if (uw == -1) {
      coefs[0] = coef_faces[0][f];
    } else {
      coefs[0] = coef_cells[0][uw];
    }

    // dw coef
    if (dw == -1) {
      coefs[1] = coef_faces[0][f];
    } else {
      coefs[1] = coef_cells[0][dw];
    }

    bool negative_coef = (coefs[0] < 0.0) || (coefs[1] < 0.0);
    AMANZI_ASSERT(!negative_coef);

    // Determine the size of the overlap region, a smooth transition region
    // near zero flux
    double flow_eps = flux_eps_;

    //Fixed coefficient in the scaling of the arithmetic mean
    double amean_order_of_supression = 15.0;

    // Determine the coefficient
    if (dw == -1)
      coef_faces[0][f] = coefs[1];
    else if (uw == -1)
      coef_faces[0][f] = coefs[0];
    else {
      double dist[2];
      dist[0] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(uw));
      dist[1] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(dw));

      double hmean = 0.0;
      if ((coefs[0] != 0.0) && (coefs[1] != 0.0))
        hmean = (dist[0] + dist[1]) / (dist[0] / coefs[0] + dist[1] / coefs[1]);

      double coef_face = hmean;
      double amean = (dist[0] * coefs[0] + dist[1] * coefs[1]) / (dist[0] + dist[1]);

      double coef_jump = 0.0;
      double amean_scaling[2];
      if (coefs[0] != coefs[1]) {
        amean_scaling[0] = (coefs[0] > 1e-15) ? std::pow(10.0,
                                                         -amean_order_of_supression * coefs[1] *
                                                           (coefs[0] + coefs[1]) /
                                                           std::pow(coefs[0] - coefs[1], 2.0)) :
                                                0.0;
        amean_scaling[1] = (coefs[1] > 1e-15) ? std::pow(10.0,
                                                         -amean_order_of_supression * coefs[0] *
                                                           (coefs[0] + coefs[1]) /
                                                           std::pow(coefs[0] - coefs[1], 2.0)) :
                                                0.0;

        coef_face += amean * amean_scaling[0];
        coef_jump = amean * std::abs(amean_scaling[0] - amean_scaling[1]);
      }

      if ((std::abs(flux_v[0][f]) < flow_eps) && (coef_jump > 1e-15)) {
        double param = std::abs(flux_v[0][f]) / flow_eps;
        double alt_coef_face = hmean + amean * amean_scaling[1];
        coef_faces[0][f] = param * coef_face + (1 - param) * alt_coef_face;
      } else
        coef_faces[0][f] = coef_face;
    }
  }
};


void
UpwindFluxHarmonicMean::UpdateDerivatives(
  const Teuchos::Ptr<State>& S,
  std::string potential_key,
  const CompositeVector& dconductivity,
  const std::vector<int>& bc_markers,
  const std::vector<double>& bc_values,
  std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>>* Jpp_faces) const
{
  AMANZI_ASSERT(0);
}
} // namespace Operators
} // namespace Amanzi
