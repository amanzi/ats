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

#include "CompositeVector.hh"
#include "State.hh"
#include "upwind_arithmetic_mean.hh"

namespace Amanzi {
namespace Operators {

UpwindArithmeticMean::UpwindArithmeticMean(const std::string& pkname, const Tag& tag)
  : pkname_(pkname), tag_(tag)
{}

void
UpwindArithmeticMean::Update(const CompositeVector& cells,
                             CompositeVector& faces,
                             const State& S,
                             const Teuchos::Ptr<Debugger>& db) const
{
  CalculateCoefficientsOnFaces(cells, "cell", faces, "face");
};

void
UpwindArithmeticMean::Update(const CompositeVector& cells,
                             const std::string cell_component,
                             CompositeVector& faces,
                             const std::string face_component,
                             const State& S,
                             const Teuchos::Ptr<Debugger>& db) const
{
  CalculateCoefficientsOnFaces(cells, cell_component, faces, face_component);
};


void
UpwindArithmeticMean::CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                                   const std::string cell_component,
                                                   CompositeVector& face_coef,
                                                   const std::string face_component) const
{
  // initialize the face coefficients
  face_coef.getComponent(face_component, true)->putScalar(0.0);
  if (face_coef.hasComponent("cell")) { face_coef.getComponent("cell", true)->putScalar(1.0); }

  // communicate ghosted cells
  cell_coef.scatterMasterToGhosted(cell_component);
  const AmanziMesh::MeshCache& m = face_coef.getMesh()->getCache();
  {
    auto face_coef_f = face_coef.viewComponent(face_component, false);
    auto cell_coef_c = cell_coef.viewComponent(cell_component, true);

    Kokkos::parallel_for(
      "upwind_arithmetic_mean", face_coef_f.extent(0), KOKKOS_LAMBDA(const int& f) {
        auto cells = m.getFaceCells(f);
        for (const auto& c : cells) { face_coef_f(f, 0) += cell_coef_c(c, 0); }
        face_coef_f(f, 0) /= cells.size();
      });
  }
  face_coef.scatterMasterToGhosted(face_component);
};


// void
// UpwindArithmeticMean::UpdateDerivatives(
//   const Teuchos::Ptr<State>& S,
//   Key potential_key,
//   const CompositeVector& dconductivity,
//   const std::vector<int>& bc_markers,
//   const std::vector<double>& bc_values,

//   std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>>* Jpp_faces) const
// {
//   // Grab derivatives
//   dconductivity.scatterMasterToGhosted("cell");
//   const Epetra_MultiVector& dcell_v = *dconductivity.viewComponent("cell", true);

//   // Grab potential
//   auto keytag = Keys::splitKeyTag(potential_key);
//   const CompositeVector& pres = S->Get<CompositeVector>(keytag.first, keytag.second);
//   pres.scatterMasterToGhosted("cell");
//   const Epetra_MultiVector& pres_v = *pres.viewComponent("cell", true);

//   // Grab mesh and allocate space
//   Teuchos::RCP<const AmanziMesh::Mesh> mesh = pres.getMesh();
//   unsigned int nfaces_owned =
//     mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
//   Jpp_faces->resize(nfaces_owned);

//   // workspace
//   double dK_dp[2];
//   double p[2];

//   for (unsigned int f = 0; f != nfaces_owned; ++f) {
//     // get neighboring cells
//     AmanziMesh::Entity_ID_List cells;
//     cells = mesh->getFaceCells(f);
//     int mcells = cells.size();

//     // create the local matrix
//     Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>> Jpp =
//       Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double>(mcells, mcells));
//     (*Jpp_faces)[f] = Jpp;

//     if (mcells == 1) {
//       if (bc_markers[f] == Operators::OPERATOR_BC_DIRICHLET) {
//         p[0] = pres_v[0][cells[0]];
//         p[1] = bc_values[f];
//         double dp = p[0] - p[1];

//         (*Jpp)(0, 0) = dp * mesh->getFaceArea(f) * dcell_v[0][cells[0]];
//       } else {
//         (*Jpp)(0, 0) = 0.;
//       }
//     } else {
//       p[0] = pres_v[0][cells[0]];
//       p[1] = pres_v[0][cells[1]];

//       dK_dp[0] = 0.5 * dcell_v[0][cells[0]];
//       dK_dp[1] = 0.5 * dcell_v[0][cells[1]];

//       (*Jpp)(0, 0) = (p[0] - p[1]) * mesh->getFaceArea(f) * dK_dp[0];
//       (*Jpp)(0, 1) = (p[0] - p[1]) * mesh->getFaceArea(f) * dK_dp[1];
//       (*Jpp)(1, 0) = -(*Jpp)(0, 0);
//       (*Jpp)(1, 1) = -(*Jpp)(0, 1);
//     }
//   }
// }


} // namespace Operators
} // namespace Amanzi
