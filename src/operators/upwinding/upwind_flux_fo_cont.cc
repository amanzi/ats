
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
#include "Epetra_IntVector.h"

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
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.Mesh();

  // initialize the face coefficients
  if (face_coef.HasComponent("cell")) { face_coef.ViewComponent("cell", true)->PutScalar(1.0); }

  // communicate needed ghost values
  cell_coef.ScatterMasterToGhosted("cell");
  slope.ScatterMasterToGhosted("cell");
  manning_coef.ScatterMasterToGhosted("cell");
  elevation.ScatterMasterToGhosted("cell");

  // pull out vectors
  const Epetra_MultiVector& flux_v = *flux.ViewComponent("face", false);
  Epetra_MultiVector& coef_faces = *face_coef.ViewComponent("face", false);
  const Epetra_MultiVector& pd_cells = *cell_coef.ViewComponent("cell", true);
  const Epetra_MultiVector& slope_v = *slope.ViewComponent("cell", false);
  const Epetra_MultiVector& manning_coef_v = *manning_coef.ViewComponent("cell", false);
  const Epetra_MultiVector& elevation_v = *elevation.ViewComponent("cell", false);
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
  double pds[2];

  int nfaces = face_coef.size("face", false);
  for (int f = 0; f != nfaces; ++f) {
    int uw = upwind_cell[f];
    int dw = downwind_cell[f];
    AMANZI_ASSERT(!((uw == -1) && (dw == -1)));

    double denominator = 0.0;
    // uw coef
    if (uw == -1) {
      denominator =
        manning_coef_v[0][dw] * std::sqrt(std::max(slope_v[0][dw], slope_regularization));
      pds[0] = coef_faces[0][f];
    } else {
      pds[0] = pd_cells[0][uw];
    }

    // dw coef
    if (dw == -1) {
      denominator =
        manning_coef_v[0][uw] * std::sqrt(std::max(slope_v[0][uw], slope_regularization));
      pds[1] = coef_faces[0][f];
    } else {
      pds[1] = pd_cells[0][dw];
    }

    if ((uw != -1) && (dw != -1)) {
      double denom[2];
      denom[0] = manning_coef_v[0][uw] * std::sqrt(std::max(slope_v[0][uw], slope_regularization));
      denom[1] = manning_coef_v[0][dw] * std::sqrt(std::max(slope_v[0][dw], slope_regularization));
      double dist[2];
      dist[0] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(uw));
      dist[1] = AmanziGeometry::norm(mesh->getFaceCentroid(f) - mesh->getCellCentroid(dw));

      denominator = (dist[0] + dist[1]) / (dist[0] / denom[0] + dist[1] / denom[1]);
    }

    double pdf = 0.0;
    // Determine the coefficient
    if (dw == -1)
      pdf = pds[1];
    else if (uw == -1)
      pdf = pds[0];
    else
      pdf = pds[0] + elevation_v[0][uw] - std::max(elevation_v[0][uw], elevation_v[0][dw]);

    double exponent = manning_exp_ + 1.0;
    coef_faces[0][f] = std::pow(std::max(pdf, 0.), exponent) / denominator;
  }
}


} // namespace Operators
} // namespace Amanzi
