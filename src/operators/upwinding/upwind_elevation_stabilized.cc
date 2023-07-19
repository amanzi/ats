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
#include "Epetra_IntVector.h"

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
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.Mesh();

  // initialize the face coefficients
  if (face_coef.HasComponent("cell")) { face_coef.ViewComponent("cell", true)->PutScalar(1.0); }

  // communicate needed ghost values
  slope.ScatterMasterToGhosted("cell");
  elevation.ScatterMasterToGhosted("cell");
  ponded_depth.ScatterMasterToGhosted("cell");
  manning_coef.ScatterMasterToGhosted("cell");
  density.ScatterMasterToGhosted("cell");

  // pull out vectors
  Epetra_MultiVector& coef_faces = *face_coef.ViewComponent("face", false);
  const Epetra_MultiVector& slope_v = *slope.ViewComponent("cell", false);
  const Epetra_MultiVector& elev_v = *elevation.ViewComponent("cell", false);
  const Epetra_MultiVector& pd_v = *ponded_depth.ViewComponent("cell", false);
  const Epetra_MultiVector& manning_coef_v = *manning_coef.ViewComponent("cell", true);
  const Epetra_MultiVector& dens_v = *density.ViewComponent("cell", true);

  const Epetra_MultiVector& elev_bf = *elevation.ViewComponent("boundary_face", false);
  const Epetra_MultiVector& pd_bf = *ponded_depth.ViewComponent("boundary_face", false);

  double slope_regularization = slope_regularization_;
  double manning_exp = manning_exp_;

  const auto& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);
  const auto& bface_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);

  // Determine the face coefficient of local faces.
  //
  // Here we upwind the h, then calculate upwinded(h+z)-max(z_uw, z_dw) as the
  // h used on the face with the usual manning model.
  //
  // Note that the upwind value here is assumed to be the max of h+z.  This is
  // always true for FV, maybe not for MFD.
  int nfaces = face_coef.size("face", false);
  for (int f = 0; f != nfaces; ++f) {
    auto fcells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    AMANZI_ASSERT(fcells.size() > 0);

    double denom[2] = { 0., 0. };
    double weight[2] = { 0., 0. };
    double pres_elev[2] = { 0., 0. };
    double elev[2] = { 0., 0. };
    double dens[2] = { 0., 0. };

    weight[0] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(fcells[0]));
    denom[0] = manning_coef_v[0][fcells[0]] *
               std::sqrt(std::max(slope_v[0][fcells[0]], slope_regularization));
    pres_elev[0] = pd_v[0][fcells[0]] + elev_v[0][fcells[0]];
    elev[0] = elev_v[0][fcells[0]];
    dens[0] = dens_v[0][fcells[0]];

    if (fcells.size() > 1) {
      weight[1] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(fcells[1]));
      denom[1] = manning_coef_v[0][fcells[1]] *
                 std::sqrt(std::max(slope_v[0][fcells[1]], slope_regularization));
      pres_elev[1] = pd_v[0][fcells[1]] + elev_v[0][fcells[1]];
      elev[1] = elev_v[0][fcells[1]];
      dens[1] = dens_v[0][fcells[1]];
    } else {
      // boundary face
      weight[1] = weight[0];
      denom[1] = denom[0];
      int bf = bface_map.LID(face_map.GID(f));
      pres_elev[1] = pd_bf[0][bf] + elev_bf[0][bf];
      elev[1] = elev_bf[0][bf];
      dens[1] = dens[0];
    }

    // harmonic mean of the denominator
    double denom_f = (weight[0] + weight[1]) / (weight[0] / denom[0] + weight[1] / denom[1]);
    double dens_f = (weight[0] + weight[1]) / (weight[0] / dens[0] + weight[1] / dens[1]);
    double h_f = std::max(pres_elev[0], pres_elev[1]) - std::max(elev[0], elev[1]);
    coef_faces[0][f] = dens_f * std::pow(h_f, 1 + manning_exp) / denom_f;
  }
};


} // namespace Operators
} // namespace Amanzi
