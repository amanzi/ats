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

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "State.hh"
#include "Debugger.hh"
#include "VerboseObject.hh"
#include "upwind_flux_split_denominator.hh"

namespace Amanzi {
namespace Operators {

UpwindFluxSplitDenominator::UpwindFluxSplitDenominator(const std::string& pkname,
                                                       const Tag& tag,
                                                       const std::string& flux,
                                                       const std::string& slope,
                                                       const std::string& manning_coef,
                                                       double flux_eps,
                                                       double slope_regularization)
  : pkname_(pkname),
    tag_(tag),
    flux_(flux),
    slope_(slope),
    manning_coef_(manning_coef),
    flux_eps_(flux_eps),
    slope_regularization_(slope_regularization)
{}


void
UpwindFluxSplitDenominator::Update(const CompositeVector& cells,
                                   CompositeVector& faces,
                                   const State& S,
                                   const Teuchos::Ptr<Debugger>& db) const
{
  Teuchos::RCP<const CompositeVector> flux = S.GetPtr<CompositeVector>(flux_, tag_);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_, tag_);
  Teuchos::RCP<const CompositeVector> manning_coef = S.GetPtr<CompositeVector>(manning_coef_, tag_);

  CalculateCoefficientsOnFaces(cells, *flux, *slope, *manning_coef, faces, db);
};


void
UpwindFluxSplitDenominator::CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                                         const CompositeVector& flux,
                                                         const CompositeVector& slope,
                                                         const CompositeVector& manning_coef,
                                                         CompositeVector& face_coef,
                                                         const Teuchos::Ptr<Debugger>& db) const
{
  // initialize the face coefficients
  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  // communicate needed ghost values
  cell_coef.scatterMasterToGhosted("cell");
  slope.scatterMasterToGhosted("cell");
  manning_coef.scatterMasterToGhosted("cell");

  // pull out vectors
  {
    const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();
    const auto flux_v = flux.viewComponent("face", false);
    auto coef_faces = face_coef.viewComponent("face", false);
    const auto coef_cells = cell_coef.viewComponent("cell", true);
    const auto slope_v = slope.viewComponent("cell", false);
    const auto manning_coef_v = manning_coef.viewComponent("cell", true);
    double slope_regularization = slope_regularization_;

    int nfaces_local = flux_v.extent(0);
    int ncells = coef_cells.extent(0);

    // Determine the face coefficient of local faces.
    // These parameters may be key to a smooth convergence rate near zero flux.
    double flow_eps = flux_eps_;

    Kokkos::parallel_for(
      "upwind_flux_split_denominator", nfaces_local, KOKKOS_LAMBDA(const int& f) {
        auto fcells = m.getFaceCells(f);

        double denominator = 0.0;
        double coefs[2] = { 0., 0. };
        double weight[2] = { 0., 0. };
        double pd[2] = { 0., 0. };
        double denom[2] = { 0., 0. };

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

        // uw coef
        if (uw == -1) {
          denominator =
            manning_coef_v(dw, 0) * Kokkos::sqrt(Kokkos::max(slope_v(dw, 0), slope_regularization));
          assert(denominator > 0);

          coefs[0] = coef_faces(f, 0) * denominator;
          coefs[1] = coef_cells(dw, 0) * denominator;

          // weighted by path length
          weight[1] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(dw));
          weight[0] = weight[1];

        } else if (dw == -1) {
          denominator =
            manning_coef_v(uw, 0) * Kokkos::sqrt(Kokkos::max(slope_v(uw, 0), slope_regularization));
          assert(denominator > 0);

          coefs[0] = coef_cells(uw, 0) * denominator;
          coefs[1] = coef_cells(uw, 0) * denominator; // downwind boundary face not defined always
          //coefs[1] = coef_faces(f,0) * denominator;

          weight[0] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(uw));
          weight[1] = weight[0];

        } else {
          assert(manning_coef_v(uw, 0) > 0);
          assert(manning_coef_v(dw, 0) > 0);
          denom[0] =
            manning_coef_v(uw, 0) * Kokkos::sqrt(Kokkos::max(slope_v(uw, 0), slope_regularization));
          denom[1] =
            manning_coef_v(dw, 0) * Kokkos::sqrt(Kokkos::max(slope_v(dw, 0), slope_regularization));

          coefs[0] = coef_cells(uw, 0) * denom[0];
          coefs[1] = coef_cells(dw, 0) * denom[1];

          // harmonic mean of the denominator
          weight[0] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(uw));
          weight[1] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(dw));
          assert(denom[0] > 0);
          assert(denom[1] > 0);
          assert(weight[0] > 0);
          assert(weight[1] > 0);
          denominator = (weight[0] + weight[1]) / (weight[0] / denom[0] + weight[1] / denom[1]);
          assert(denominator > 0);
        }

        // Determine the coefficient
        assert(denominator > 0);
        assert(coefs[0] >= 0 && coefs[1] >= 0);
        if (coefs[1] > coefs[0]) {
          // downwind ponded depth is larger
          if ((coefs[0] != 0.0)) {
            // harmonic mean (smoothly approaches zero for upwind coef = 0)
            coef_faces(f, 0) =
              (weight[0] + weight[1]) / (weight[0] / coefs[0] + weight[1] / coefs[1]);
          } else {
            // harmonic mean of zero is zero
            coef_faces(f, 0) = 0.0;
          }
        } else if (Kokkos::abs(flux_v(f, 0)) >= flow_eps) {
          // upwind ponded depth is larger, flux potential is nonzero
          // arithmetic mean (smoothly stays nonzero as downwind coef approches zero)
          coef_faces(f, 0) =
            (weight[0] * coefs[0] + weight[1] * coefs[1]) / (weight[0] + weight[1]);
        } else {
          // upwind ponded depth is larger, flux potential approaches zero
          // smoothly vary between harmonic and arithmetic means
          double param = Kokkos::abs(flux_v(f, 0)) / flow_eps;
          double amean = (weight[0] * coefs[0] + weight[1] * coefs[1]) / (weight[0] + weight[1]);
          double hmean = 0.0;
          if ((coefs[0] != 0.0) && (coefs[1] != 0.0))
            hmean = (weight[0] + weight[1]) / (weight[0] / coefs[0] + weight[1] / coefs[1]);
          coef_faces(f, 0) = param * amean + (1 - param) * hmean;
        }

        // divide by harmonic mean denominator
        coef_faces(f, 0) /= denominator;
      });
  }
  face_coef.scatterMasterToGhosted("face");
}

} // namespace Operators
} // namespace Amanzi
