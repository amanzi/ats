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
#include "Epetra_IntVector.h"

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
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.Mesh();

  // initialize the face coefficients
  if (face_coef.HasComponent("cell")) { face_coef.ViewComponent("cell", true)->PutScalar(1.0); }

  // communicate needed ghost values
  cell_coef.ScatterMasterToGhosted("cell");
  slope.ScatterMasterToGhosted("cell");
  manning_coef.ScatterMasterToGhosted("cell");

  // pull out vectors
  const Epetra_MultiVector& flux_v = *flux.ViewComponent("face", false);
  Epetra_MultiVector& coef_faces = *face_coef.ViewComponent("face", false);
  const Epetra_MultiVector& coef_cells = *cell_coef.ViewComponent("cell", true);
  const Epetra_MultiVector& slope_v = *slope.ViewComponent("cell", false);
  const Epetra_MultiVector& manning_coef_v = *manning_coef.ViewComponent("cell", true);
  double slope_regularization = slope_regularization_;

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
  //  double flow_eps_factor = 1.;
  //  double min_flow_eps = 1.e-8;
  int nfaces = face_coef.size("face", false);
  for (int f = 0; f != nfaces; ++f) {
    int uw = upwind_cell[f];
    int dw = downwind_cell[f];
    AMANZI_ASSERT(!((uw == -1) && (dw == -1)));

    double denominator = 0.0;
    double coefs[2] = { 0., 0. };
    double weight[2] = { 0., 0. };
    double pd[2] = { 0., 0. };
    double denom[2] = { 0., 0. };

    // uw coef
    if (uw == -1) {
      denominator =
        manning_coef_v[0][dw] * std::sqrt(std::max(slope_v[0][dw], slope_regularization));
      AMANZI_ASSERT(denominator > 0);

      coefs[0] = coef_faces[0][f] * denominator;
      coefs[1] = coef_cells[0][dw] * denominator;

      // weighted by path length
      weight[1] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(dw));
      weight[0] = weight[1];

    } else if (dw == -1) {
      denominator =
        manning_coef_v[0][uw] * std::sqrt(std::max(slope_v[0][uw], slope_regularization));
      AMANZI_ASSERT(denominator > 0);

      coefs[0] = coef_cells[0][uw] * denominator;
      coefs[1] = coef_cells[0][uw] * denominator; // downwind boundary face not defined always
      //coefs[1] = coef_faces[0][f] * denominator;

      weight[0] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(uw));
      weight[1] = weight[0];

    } else {
      AMANZI_ASSERT(manning_coef_v[0][uw] > 0);
      AMANZI_ASSERT(manning_coef_v[0][dw] > 0);
      denom[0] = manning_coef_v[0][uw] * std::sqrt(std::max(slope_v[0][uw], slope_regularization));
      denom[1] = manning_coef_v[0][dw] * std::sqrt(std::max(slope_v[0][dw], slope_regularization));

      coefs[0] = coef_cells[0][uw] * denom[0];
      coefs[1] = coef_cells[0][dw] * denom[1];

      // harmonic mean of the denominator
      weight[0] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(uw));
      weight[1] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(dw));
      AMANZI_ASSERT(denom[0] > 0);
      AMANZI_ASSERT(denom[1] > 0);
      AMANZI_ASSERT(weight[0] > 0);
      AMANZI_ASSERT(weight[1] > 0);
      denominator = (weight[0] + weight[1]) / (weight[0] / denom[0] + weight[1] / denom[1]);
      AMANZI_ASSERT(denominator > 0);
    }

    double flow_eps = flux_eps_;

    // Determine the coefficient
    AMANZI_ASSERT(denominator > 0);
    AMANZI_ASSERT(coefs[0] >= 0 && coefs[1] >= 0);
    if (coefs[1] > coefs[0]) {
      // downwind ponded depth is larger
      if ((coefs[0] != 0.0)) {
        // harmonic mean (smoothly approaches zero for upwind coef = 0)
        coef_faces[0][f] = (weight[0] + weight[1]) / (weight[0] / coefs[0] + weight[1] / coefs[1]);
      } else {
        // harmonic mean of zero is zero
        coef_faces[0][f] = 0.0;
      }
    } else if (std::abs(flux_v[0][f]) >= flow_eps) {
      // upwind ponded depth is larger, flux potential is nonzero
      // arithmetic mean (smoothly stays nonzero as downwind coef approches zero)
      coef_faces[0][f] = (weight[0] * coefs[0] + weight[1] * coefs[1]) / (weight[0] + weight[1]);
    } else {
      // upwind ponded depth is larger, flux potential approaches zero
      // smoothly vary between harmonic and arithmetic means
      double param = std::abs(flux_v[0][f]) / flow_eps;
      double amean = (weight[0] * coefs[0] + weight[1] * coefs[1]) / (weight[0] + weight[1]);
      double hmean = 0.0;
      if ((coefs[0] != 0.0) && (coefs[1] != 0.0))
        hmean = (weight[0] + weight[1]) / (weight[0] / coefs[0] + weight[1] / coefs[1]);
      coef_faces[0][f] = param * amean + (1 - param) * hmean;
    }

    // divide by harmonic mean denominator
    coef_faces[0][f] /= denominator;
  }
};


} // namespace Operators
} // namespace Amanzi
