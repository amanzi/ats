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
#include "upwind_elevation_stabilized.hh"

namespace Amanzi {
namespace Operators {

UpwindElevationStabilized::UpwindElevationStabilized(const std::string& pkname,
                                                     const Tag& tag,
                                                     const std::string& slope,
                                                     const std::string& manning_coef,
                                                     const std::string& ponded_depth,
                                                     const std::string& elevation,
                                                     const std::string& density,
                                                     double slope_regularization,
                                                     double manning_exp)
  : pkname_(pkname),
    tag_(tag),
    slope_(slope),
    manning_coef_(manning_coef),
    ponded_depth_(ponded_depth),
    elevation_(elevation),
    density_(density),
    slope_regularization_(slope_regularization),
    manning_exp_(manning_exp)
{}


void
UpwindElevationStabilized::Update(const CompositeVector& cells,
                                  CompositeVector& faces,
                                  const State& S,
                                  const Teuchos::Ptr<Debugger>& db) const
{
  Teuchos::RCP<const CompositeVector> slope = S.GetPtr<CompositeVector>(slope_, tag_);
  Teuchos::RCP<const CompositeVector> elev = S.GetPtr<CompositeVector>(elevation_, tag_);
  Teuchos::RCP<const CompositeVector> pd = S.GetPtr<CompositeVector>(ponded_depth_, tag_);
  Teuchos::RCP<const CompositeVector> manning_coef = S.GetPtr<CompositeVector>(manning_coef_, tag_);
  Teuchos::RCP<const CompositeVector> density = S.GetPtr<CompositeVector>(density_, tag_);

  CalculateCoefficientsOnFaces(*slope, *manning_coef, *pd, *elev, *density, faces, db);
};


void
UpwindElevationStabilized::CalculateCoefficientsOnFaces(const CompositeVector& slope,
                                                        const CompositeVector& manning_coef,
                                                        const CompositeVector& ponded_depth,
                                                        const CompositeVector& elevation,
                                                        const CompositeVector& density,
                                                        CompositeVector& face_coef,
                                                        const Teuchos::Ptr<Debugger>& db) const
{
  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();

  // initialize the face coefficients
  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  // communicate needed ghost values
  slope.scatterMasterToGhosted("cell");
  elevation.scatterMasterToGhosted("cell");
  ponded_depth.scatterMasterToGhosted("cell");
  manning_coef.scatterMasterToGhosted("cell");
  density.scatterMasterToGhosted("cell");

  {
    // pull out vectors
    auto coef_faces = face_coef.viewComponent("face", false);
    const auto slope_v = slope.viewComponent("cell", false);
    const auto elev_v = elevation.viewComponent("cell", false);
    const auto pd_v = ponded_depth.viewComponent("cell", false);
    const auto manning_coef_v = manning_coef.viewComponent("cell", true);
    const auto dens_v = density.viewComponent("cell", true);

    const auto elev_bf = elevation.viewComponent("boundary_face", false);
    const auto pd_bf = ponded_depth.viewComponent("boundary_face", false);

    double slope_regularization = slope_regularization_;
    double manning_exp = manning_exp_ + 1;

    // const auto& face_map = m.getMap(AmanziMesh::Entity_kind::FACE,false);
    // const auto& bface_map = m.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);

    // Determine the face coefficient of local faces.
    //
    // Here we upwind the h, then calculate upwinded(h+z)-max(z_uw, z_dw) as the
    // h used on the face with the usual manning model.
    //
    // Note that the upwind value here is assumed to be the max of h+z.  This is
    // always true for FV, maybe not for MFD.
    int nfaces_local = coef_faces.extent(0);
    Kokkos::parallel_for(
      "upwind_flux_elevation_stabilized", nfaces_local, KOKKOS_LAMBDA(const int& f) {
        auto fcells = m.getFaceCells(f);

        double denom[2] = { 0., 0. };
        double weight[2] = { 0., 0. };
        double pres_elev[2] = { 0., 0. };
        double elev[2] = { 0., 0. };
        double dens[2] = { 0., 0. };

        weight[0] =
          AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(fcells[0]));
        denom[0] = manning_coef_v(fcells[0], 0) *
                   std::sqrt(fmax(slope_v(fcells[0], 0), slope_regularization));
        pres_elev[0] = pd_v(fcells[0], 0) + elev_v(fcells[0], 0);
        elev[0] = elev_v(fcells[0], 0);
        dens[0] = dens_v(fcells[0], 0);

        if (fcells.size() > 1) {
          weight[1] =
            AmanziGeometry::norm(m.getFaceCentroid(f) - m.getCellCentroid(fcells[1]));
          denom[1] = manning_coef_v(fcells[1], 0) *
                     std::sqrt(fmax(slope_v(fcells[1], 0), slope_regularization));
          pres_elev[1] = pd_v(fcells[1], 0) + elev_v(fcells[1], 0);
          elev[1] = elev_v(fcells[1], 0);
          dens[1] = dens_v(fcells[1], 0);
        } else {
          // boundary face
          weight[1] = weight[0];
          denom[1] = denom[0];

          // NOTE: THIS CURRENTLY BREAKS THINGS>>> not sure how to deal with this
          int bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(m, f);
          // END NOTE
          pres_elev[1] = pd_bf(bf, 0) + elev_bf(bf, 0);
          elev[1] = elev_bf(bf, 0);
          dens[1] = dens[0];
        }

        // harmonic mean of the denominator
        double denom_f = (weight[0] + weight[1]) / (weight[0] / denom[0] + weight[1] / denom[1]);
        double dens_f = (weight[0] + weight[1]) / (weight[0] / dens[0] + weight[1] / dens[1]);
        double h_f = fmax(pres_elev[0], pres_elev[1]) - fmax(elev[0], elev[1]);
        coef_faces(f, 0) = dens_f * pow(h_f, manning_exp) / denom_f;
      });
  }
  face_coef.scatterMasterToGhosted("face");
};


} // namespace Operators
} // namespace Amanzi
