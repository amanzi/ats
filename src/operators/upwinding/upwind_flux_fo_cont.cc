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
#include "upwind_flux_fo_cont.hh"

namespace Amanzi {
namespace Operators {

UpwindFluxFOCont::UpwindFluxFOCont(const std::string& pkname,
                                   const Tag& tag,
                                   const Key& flux,
                                   const Key& slope,
                                   const Key& manning_coef,
                                   const Key& elevation,
                                   double slope_regularization,
                                   double manning_exp)
  : pkname_(pkname),
    tag_(tag),
    flux_(flux),
    slope_(slope),
    manning_coef_(manning_coef),
    elevation_(elevation),
    slope_regularization_(slope_regularization),
    manning_exp_(manning_exp){};


void
UpwindFluxFOCont::Update(const CompositeVector& cells,
                         CompositeVector& faces,
                         const State& S,
                         const Teuchos::Ptr<Debugger>& db) const
{
  Teuchos::RCP<const CompositeVector> flux = S.GetPtr<CompositeVector>(flux_, tag_);
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_, tag_);
  Teuchos::RCP<const CompositeVector> manning_coef = S.GetPtr<CompositeVector>(manning_coef_, tag_);
  Teuchos::RCP<const CompositeVector> elevation = S.GetPtr<CompositeVector>(elevation_, tag_);
  CalculateCoefficientsOnFaces(cells, *flux, *slope, *manning_coef, *elevation, faces, db);
};


void
UpwindFluxFOCont::CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                               const CompositeVector& flux,
                                               const CompositeVector& slope,
                                               const CompositeVector& manning_coef,
                                               const CompositeVector& elevation,
                                               CompositeVector& face_coef,
                                               const Teuchos::Ptr<Debugger>& db) const
{
  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();

  // initialize the face coefficients
  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  // communicate needed ghost values
  cell_coef.scatterMasterToGhosted("cell");
  slope.scatterMasterToGhosted("cell");
  manning_coef.scatterMasterToGhosted("cell");
  elevation.scatterMasterToGhosted("cell");

  {
    // pull out vectors
    const auto flux_v = flux.viewComponent("face", false);
    auto coef_faces = face_coef.viewComponent("face", false);
    const auto pd_cells = cell_coef.viewComponent("cell", true);
    const auto slope_v = slope.viewComponent("cell", false);
    const auto manning_coef_v = manning_coef.viewComponent("cell", false);
    const auto elevation_v = elevation.viewComponent("cell", false);
    double slope_regularization = slope_regularization_;

    int nfaces_local = flux_v.extent(0);

    // Determine the face coefficient of local faces.
    // These parameters may be key to a smooth convergence rate near zero flux.
    double exponent = manning_exp_ + 1.0;
    Kokkos::parallel_for(
      "upwind_flux_fo_cont", nfaces_local, KOKKOS_LAMBDA(const int& f) {
        auto fcells = m.getFaceCells(f);

        double pds[2] = { 0., 0. };

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

        double denominator = 0.0;
        // uw coef
        if (uw == -1) {
          denominator =
            manning_coef_v(dw, 0) * sqrt(fmax(slope_v(dw, 0), slope_regularization));
          pds[0] = coef_faces(f, 0);
        } else {
          pds[0] = pd_cells(uw, 0);
        }

        // dw coef
        if (dw == -1) {
          denominator =
            manning_coef_v(uw, 0) * sqrt(fmax(slope_v(uw, 0), slope_regularization));
          pds[1] = coef_faces(f, 0);
        } else {
          pds[1] = pd_cells(dw, 0);
        }

        if ((uw != -1) && (dw != -1)) {
          double denom[2];
          denom[0] =
            manning_coef_v(uw, 0) * sqrt(fmax(slope_v(uw, 0), slope_regularization));
          denom[1] =
            manning_coef_v(dw, 0) * sqrt(fmax(slope_v(dw, 0), slope_regularization));
          double dist[2];
          dist[0] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(uw));
          dist[1] = AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(dw));
          denominator = (dist[0] + dist[1]) / (dist[0] / denom[0] + dist[1] / denom[1]);
        }

        double pdf = 0.0;
        // Determine the coefficient
        if (dw == -1)
          pdf = pds[1];
        else if (uw == -1)
          pdf = pds[0];
        else
          pdf = pds[0] + elevation_v(uw, 0) - fmax(elevation_v(uw, 0), elevation_v(dw, 0));

        coef_faces(f, 0) = pow(fmax(pdf, 0.), exponent) / denominator;
      });
  }
  face_coef.scatterMasterToGhosted("face");
}

} // namespace Operators
} // namespace Amanzi
