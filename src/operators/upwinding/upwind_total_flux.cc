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
#include "upwind_total_flux.hh"

namespace Amanzi {
namespace Operators {

UpwindTotalFlux::UpwindTotalFlux(const std::string& pkname,
                                 const Tag& tag,
                                 const std::string& flux,
                                 double flux_eps)
  : pkname_(pkname), tag_(tag), flux_(flux), flux_eps_(flux_eps){};


void
UpwindTotalFlux::Update(const CompositeVector& cells,
                        CompositeVector& faces,
                        const State& S,
                        const Teuchos::Ptr<Debugger>& db) const
{
  Teuchos::RCP<const CompositeVector> flux = S.GetPtr<CompositeVector>(flux_, tag_);
  CalculateCoefficientsOnFaces(cells, "cell", *flux, faces, "face", db);
};

void
UpwindTotalFlux::Update(const CompositeVector& cells,
                        const std::string cell_component,
                        CompositeVector& faces,
                        const std::string face_component,
                        const State& S,
                        const Teuchos::Ptr<Debugger>& db) const
{
  Teuchos::RCP<const CompositeVector> flux = S.GetPtr<CompositeVector>(flux_, tag_);
  CalculateCoefficientsOnFaces(cells, cell_component, *flux, faces, face_component, db);
};


void
UpwindTotalFlux::CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                              const std::string cell_component,
                                              const CompositeVector& flux,
                                              CompositeVector& face_coef,
                                              const std::string face_component,
                                              const Teuchos::Ptr<Debugger>& db) const
{
  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();

  // initialize the face coefficients
  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  // communicate needed ghost values
  cell_coef.scatterMasterToGhosted(cell_component);

  {
    // pull out vectors
    const auto flux_v = flux.viewComponent("face", false);
    auto coef_faces = face_coef.viewComponent(face_component, false);
    const auto coef_cells = cell_coef.viewComponent(cell_component, true);

    int nfaces_local = coef_faces.extent(0);
    bool has_cells = face_coef.hasComponent("cell");
    // NOTE: this capability exists in master, but not sure why!
    // CompositeVector::cView_type face_cell_coef;
    // if (has_cells) face_cell_coef = face_coef.viewComponent("cell", true);
    double flow_eps = flux_eps_;

    Kokkos::parallel_for(
      "upwind_total_flux", nfaces_local, KOKKOS_LAMBDA(const int& f) {
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

        // Determine the face coefficient of local faces.
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

        // Determine the coefficient
        if (Kokkos::abs(flux_v(f, 0)) >= flow_eps) {
          coef_faces(f, 0) = coefs[0];
        } else {
          // Parameterization of a linear scaling between upwind and downwind.
          double param = Kokkos::abs(flux_v(f, 0)) / (2 * flow_eps) + 0.5;
          // if (!(param >= 0.5) || !(param <= 1.0)) {
          //   std::cout << "BAD FLUX! on face " << f << std::endl;
          //   std::cout << "  flux = " << flux_v(f, 0) << std::endl;
          //   std::cout << "  param = " << param << std::endl;
          //   std::cout << "  flow_eps = " << flow_eps << std::endl;
          // }

          assert(param >= 0.5);
          assert(param <= 1.0);

          coef_faces(f, 0) = coefs[0] * param + coefs[1] * (1. - param);
        }
      });
  }
  face_coef.scatterMasterToGhosted("face");
};


// void
// UpwindTotalFlux::UpdateDerivatives(
//   const Teuchos::Ptr<State>& S,
//   std::string potential_key,
//   const CompositeVector& dconductivity,
//   const std::vector<int>& bc_markers,
//   const std::vector<double>& bc_values,
//   std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>>* Jpp_faces) const
// {
//   // Grab derivatives
//   dconductivity.scatterMasterToGhosted("cell");
//   const Epetra_MultiVector& dcell_v = *dconductivity.viewComponent("cell", true);

//   // Grab potential
//   Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(potential_key, tag_);
//   pres->scatterMasterToGhosted("cell");
//   const Epetra_MultiVector& pres_v = *pres->viewComponent("cell", true);

//   // Grab flux direction
//   const Epetra_MultiVector& flux_v =
//     *S->Get<CompositeVector>(flux_, tag_).viewComponent("face", false);

//   // Grab mesh and allocate space
//   Teuchos::RCP<const AmanziMesh::Mesh> mesh = dconductivity.getMesh();
//   unsigned int nfaces_owned =
//     mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
//   Jpp_faces->resize(nfaces_owned);

//   // workspace
//   double dK_dp[2];
//   double p[2];


//   // Identify upwind/downwind cells for each local face.  Note upwind/downwind
//   // may be a ghost cell.
//   Epetra_IntVector upwind_cell(mesh->getMap(AmanziMesh::Entity_kind::FACE,true));
//   upwind_cell.PutValue(-1);
//   Epetra_IntVector downwind_cell(mesh->getMap(AmanziMesh::Entity_kind::FACE,true));
//   downwind_cell.PutValue(-1);

//   AmanziMesh::Entity_ID_List faces;
//   std::vector<int> fdirs;

//   int ncells = dcell_v.MyLength();
//   for (int c = 0; c != ncells; ++c) {
//     mesh->getCellFacesAndDirections(c, &faces, &fdirs);

//     for (unsigned int n = 0; n != faces.size(); ++n) {
//       int f = faces[n];

//       if (f < nfaces_owned) {
//         if (flux_v[0][f] * fdirs[n] > 0) {
//           upwind_cell[f] = c;
//         } else if (flux_v[0][f] * fdirs[n] < 0) {
//           downwind_cell[f] = c;
//         } else {
//           // We don't care, but we have to get one into upwind and the other
//           // into downwind.
//           if (upwind_cell[f] == -1) {
//             upwind_cell[f] = c;
//           } else {
//             downwind_cell[f] = c;
//           }
//         }
//       }
//     }
//   }


//   for (unsigned int f = 0; f != nfaces_owned; ++f) {
//     int uw = upwind_cell[f];
//     int dw = downwind_cell[f];
//     AMANZI_ASSERT(!((uw == -1) && (dw == -1)));

//     AmanziMesh::Entity_ID_List cells;
//     cells = mesh->getFaceCells(f);
//     int mcells = cells.size();

//     // uw coef
//     if (uw == -1) {
//       // boundary, upwind is the boundary
//       if (std::abs(flux_v[0][f]) >= flux_eps_) {
//         // flux coming from boundary, derivs are zero
//         dK_dp[0] = 0.;
//       } else {
//         // Parameterization of a linear scaling between upwind and downwind.
//         double param = std::abs(flux_v[0][f]) / (2 * flux_eps_) + 0.5;

//         // ignoring dparam_dp... not sure how we would include that
//         dK_dp[0] = (1. - param) * dcell_v[0][dw];
//       }

//     } else if (dw == -1) {
//       // boundary, upwind is the cell
//       if (std::abs(flux_v[0][f]) >= flux_eps_) {
//         dK_dp[0] = dcell_v[0][uw];
//       } else {
//         double param = std::abs(flux_v[0][f]) / (2 * flux_eps_) + 0.5;
//         dK_dp[0] = param * dcell_v[0][uw];
//       }

//     } else {
//       // non-boundary
//       if (std::abs(flux_v[0][f]) >= flux_eps_) {
//         if (uw == cells[0]) {
//           dK_dp[0] = dcell_v[0][uw];
//           dK_dp[1] = 0.;
//         } else {
//           dK_dp[1] = dcell_v[0][uw];
//           dK_dp[0] = 0.;
//         }
//       } else {
//         double param = std::abs(flux_v[0][f]) / (2 * flux_eps_) + 0.5;
//         if (uw == cells[0]) {
//           dK_dp[0] = param * dcell_v[0][uw];
//           dK_dp[1] = (1 - param) * dcell_v[0][dw];
//         } else {
//           dK_dp[1] = param * dcell_v[0][uw];
//           dK_dp[0] = (1 - param) * dcell_v[0][dw];
//         }
//       }
//     }

//     // create the local matrix
//     Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>> Jpp =
//       Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double>(mcells, mcells));
//     (*Jpp_faces)[f] = Jpp;

//     if (mcells == 1) {
//       if (bc_markers[f] == Operators::OPERATOR_BC_DIRICHLET) {
//         // determine flux
//         p[0] = pres_v[0][cells[0]];
//         p[1] = bc_values[f];
//         double dp = p[0] - p[1];

//         (*Jpp)(0, 0) = dp * mesh->getFaceArea(f) * dK_dp[0];
//       } else {
//         (*Jpp)(0, 0) = 0.;
//       }

//     } else {
//       p[0] = pres_v[0][cells[0]];
//       p[1] = pres_v[0][cells[1]];

//       (*Jpp)(0, 0) = (p[0] - p[1]) * mesh->getFaceArea(f) * dK_dp[0];
//       (*Jpp)(0, 1) = (p[0] - p[1]) * mesh->getFaceArea(f) * dK_dp[1];
//       (*Jpp)(1, 0) = -(*Jpp)(0, 0);
//       (*Jpp)(1, 1) = -(*Jpp)(0, 1);
//     }
//   }
// }
} // namespace Operators
} // namespace Amanzi
