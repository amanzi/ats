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
  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();

  // initialize the face coefficients
  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  // communicate needed ghost values
  cell_coef.scatterMasterToGhosted("cell");

  // pull out vectors
  {
    const auto flux_v = flux.viewComponent("face", false);
    auto coef_faces = face_coef.viewComponent("face", false);
    const auto coef_cells = cell_coef.viewComponent("cell", true);
    double flow_eps = flux_eps_;

    int nfaces_local = coef_faces.extent(0);
    Kokkos::parallel_for(
      "upwind_flux_harmonic_mean", nfaces_local, KOKKOS_LAMBDA(const int& f) {
        auto fcells = m.getFaceCells(f);

        int uw = -1, dw = -1;
        int c0 = fcells(0);
        int orientation = 0;
        m.getFaceNormal(f, c0, &orientation);
        if (flux_v(f, 0) * orientation > 0) {
          uw = c0;
          if (fcells.size() == 2) dw = fcells(1);
        } else {
          dw = c0;
          if (fcells.size() == 2) uw = fcells(1);
        }
        assert(!((uw == -1) && (dw == -1)));

        double coefs[2];
        // uw coef
        if (uw == -1) {
          coefs[0] = coef_faces(f, 0);
        } else {
          coefs[0] = coef_cells(uw, 0);
        }

        // dw coef
        if (dw == -1) {
          coefs[1] = coef_faces(f, 0);
        } else {
          coefs[1] = coef_cells(dw, 0);
        }

        assert(!(coefs[0] < 0.0) || (coefs[1] < 0.0));

        // Determine the size of the overlap region, a smooth transition region
        // near zero flux

        // Fixed coefficient in the scaling of the arithmetic mean
        double amean_order_of_supression = 15.0;

        // Determine the coefficient
        if (dw == -1) {
          coef_faces(f, 0) = coefs[1];
        } else if (uw == -1) {
          coef_faces(f, 0) = coefs[0];
        } else {
          double dist[2];
          dist[0] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(uw));
          dist[1] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(dw));

          double hmean = 0.0;
          if ((coefs[0] != 0.0) && (coefs[1] != 0.0))
            hmean = (dist[0] + dist[1]) / (dist[0] / coefs[0] + dist[1] / coefs[1]);

          double coef_face = hmean;
          double amean = (dist[0] * coefs[0] + dist[1] * coefs[1]) / (dist[0] + dist[1]);

          double coef_jump = 0.0;
          double amean_scaling[2];
          if (coefs[0] != coefs[1]) {
            amean_scaling[0] = (coefs[0] > 1e-15) ? pow(10.0,
                                                             -amean_order_of_supression * coefs[1] *
                                                               (coefs[0] + coefs[1]) /
                                                               pow(coefs[0] - coefs[1], 2.0)) :
                                                    0.0;
            amean_scaling[1] = (coefs[1] > 1e-15) ? pow(10.0,
                                                             -amean_order_of_supression * coefs[0] *
                                                               (coefs[0] + coefs[1]) /
                                                               pow(coefs[0] - coefs[1], 2.0)) :
                                                    0.0;

            coef_face += amean * amean_scaling[0];
            coef_jump = amean * fabs(amean_scaling[0] - amean_scaling[1]);
          }

          if ((fabs(flux_v(f, 0)) < flow_eps) && (coef_jump > 1e-15)) {
            double param = fabs(flux_v(f, 0)) / flow_eps;
            double alt_coef_face = hmean + amean * amean_scaling[1];
            coef_faces(f, 0) = param * coef_face + (1 - param) * alt_coef_face;
          } else {
            coef_faces(f, 0) = coef_face;
          }
        }
      });
  }
  face_coef.scatterMasterToGhosted("face");
};


// void
// UpwindFluxHarmonicMean::UpdateDerivatives(
//   const Teuchos::Ptr<State>& S,
//   std::string potential_key,
//   const CompositeVector& dconductivity,
//   const std::vector<int>& bc_markers,
//   const std::vector<double>& bc_values,
//   std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>>* Jpp_faces) const
// {
//   AMANZI_ASSERT(0);
// }
} // namespace Operators
} // namespace Amanzi
