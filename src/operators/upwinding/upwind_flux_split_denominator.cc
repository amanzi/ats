/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
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

UpwindFluxSplitDenominator::UpwindFluxSplitDenominator(std::string pkname,
                                                       std::string cell_coef,
                                                       std::string face_coef,
                                                       std::string flux,
                                                       double flux_eps,
                                                       std::string slope,
                                                       std::string manning_coef,
                                                       double slope_regularization,
                                                       std::string ponded_depth) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef),
    flux_(flux),
    flux_eps_(flux_eps),
    slope_(slope),
    manning_coef_(manning_coef),
    slope_regularization_(slope_regularization),
    ponded_depth_(ponded_depth) {};


void UpwindFluxSplitDenominator::Update(const Teuchos::Ptr<State>& S,
                                        const Teuchos::Ptr<Debugger>& db) {

  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_coef_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_);
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_coef_, pkname_);

  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_);
  Teuchos::RCP<const CompositeVector> manning_coef = S->GetFieldData(manning_coef_);
  Teuchos::RCP<const CompositeVector> ponded_depth = S->GetFieldData(ponded_depth_);

  CalculateCoefficientsOnFaces(*cell, *flux, *slope, *manning_coef, *ponded_depth, face.ptr(), db);
};


void UpwindFluxSplitDenominator::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const CompositeVector& flux,
        const CompositeVector& slope,
        const CompositeVector& manning_coef,
        const CompositeVector& ponded_depth,
        const Teuchos::Ptr<CompositeVector>& face_coef,
        const Teuchos::Ptr<Debugger>& db) {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef->Mesh();

  // initialize the face coefficients
  if (face_coef->HasComponent("cell")) {
    face_coef->ViewComponent("cell",true)->PutScalar(1.0);
  }

  // communicate needed ghost values
  cell_coef.ScatterMasterToGhosted("cell");
  slope.ScatterMasterToGhosted("cell");
  manning_coef.ScatterMasterToGhosted("cell");
  ponded_depth.ScatterMasterToGhosted("cell");

  // pull out vectors
  const Epetra_MultiVector& flux_v = *flux.ViewComponent("face",false);
  Epetra_MultiVector& coef_faces = *face_coef->ViewComponent("face",false);
  const Epetra_MultiVector& coef_cells = *cell_coef.ViewComponent("cell",true);
  const Epetra_MultiVector& slope_v = *slope.ViewComponent("cell",false);
  const Epetra_MultiVector& manning_coef_v = *manning_coef.ViewComponent("cell",false);
  const Epetra_MultiVector& ponded_depth_v = *ponded_depth.ViewComponent("cell",false);
  double slope_regularization = slope_regularization_;

  // Identify upwind/downwind cells for each local face.  Note upwind/downwind
  // may be a ghost cell.
  Epetra_IntVector upwind_cell(*face_coef->ComponentMap("face",true));
  upwind_cell.PutValue(-1);
  Epetra_IntVector downwind_cell(*face_coef->ComponentMap("face",true));
  downwind_cell.PutValue(-1);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  int nfaces_local = flux.size("face",false);

  int ncells = cell_coef.size("cell",true);
  for (int c=0; c!=ncells; ++c) {
    mesh->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (unsigned int n=0; n!=faces.size(); ++n) {
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
  double coefs[2];

  int nfaces = face_coef->size("face",false);
  for (int f=0; f!=nfaces; ++f) {
    int uw = upwind_cell[f];
    int dw = downwind_cell[f];
    AMANZI_ASSERT(!((uw == -1) && (dw == -1)));

    double denominator = 0.0;
    // uw coef
    if (uw == -1) {
      denominator = manning_coef_v[0][dw] * std::sqrt(std::max(slope_v[0][dw], slope_regularization));
      coefs[0] = coef_faces[0][f];
    } else {
      coefs[0] = coef_cells[0][uw];
    }

    // dw coef
    if (dw == -1) {
      denominator = manning_coef_v[0][uw] * std::sqrt(std::max(slope_v[0][uw], slope_regularization));
      coefs[1] = coef_faces[0][f];
    } else {
      coefs[1] = coef_cells[0][dw];
    }

    if ((uw != -1)&&(dw != -1)) {
      double denom[2];
      denom[0] = manning_coef_v[0][uw] * std::sqrt(std::max(slope_v[0][uw], slope_regularization));
      denom[1] = manning_coef_v[0][dw] * std::sqrt(std::max(slope_v[0][dw], slope_regularization));
      double dist[2];
      dist[0] = AmanziGeometry::norm(mesh->face_centroid(f) - mesh->cell_centroid(uw));
      dist[1] = AmanziGeometry::norm(mesh->face_centroid(f) - mesh->cell_centroid(dw));

      denominator = (dist[0] + dist[1])/(dist[0]/denom[0] + dist[1]/denom[1]);
    }

    double flow_eps = flux_eps_;

    // Determine the coefficient
    if (dw == -1) coef_faces[0][f] = coefs[1];
    else if (uw == -1) coef_faces[0][f] = coefs[0];
    else {
      double dist[2];
      dist[0] = AmanziGeometry::norm(mesh->face_centroid(f) - mesh->cell_centroid(uw));
      dist[1] = AmanziGeometry::norm(mesh->face_centroid(f) - mesh->cell_centroid(dw));

      if (ponded_depth_v[0][dw] > ponded_depth_v[0][uw]) {
        if ((coefs[0] != 0.0) && (coefs[1] != 0.0))
          coef_faces[0][f] = (dist[0] + dist[1])/(dist[0]/coefs[0] + dist[1]/coefs[1]);
        else coef_faces[0][f] = 0.0;
      } else if (std::abs(flux_v[0][f]) >= flow_eps) {
        coef_faces[0][f] = (dist[0]*coefs[0] + dist[1]*coefs[1])/(dist[0] + dist[1]);
      }
      else {
        double param = std::abs(flux_v[0][f]) / flow_eps;
        double amean = (dist[0]*coefs[0] + dist[1]*coefs[1])/(dist[0] + dist[1]);
        double hmean = 0.0;
        if ((coefs[0] != 0.0) && (coefs[1] != 0.0))
          hmean = (dist[0] + dist[1])/(dist[0]/coefs[0] + dist[1]/coefs[1]);
        coef_faces[0][f] = param*amean + (1 - param)*hmean;
      }
    }
    coef_faces[0][f] /= denominator;
  }
};


void
UpwindFluxSplitDenominator::UpdateDerivatives(const Teuchos::Ptr<State>& S,
                                              std::string potential_key,
                                              const CompositeVector& dconductivity,
                                              const std::vector<int>& bc_markers,
                                              const std::vector<double>& bc_values,
                                              std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > >* Jpp_faces) const {
  AMANZI_ASSERT(0);
}
} //namespace
} //namespace
