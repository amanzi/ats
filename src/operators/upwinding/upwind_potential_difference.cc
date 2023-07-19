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
// Scheme for taking coefficients for div-grad operators from cells to faces.
// Upwinds based upon a potential vector, with an overlap region size
// determined by an (optional) other field.
// -----------------------------------------------------------------------------

#include "CompositeVector.hh"
#include "State.hh"
#include "upwind_potential_difference.hh"

namespace Amanzi {
namespace Operators {


UpwindPotentialDifference::UpwindPotentialDifference(const std::string& pkname,
                                                     const Tag& tag,
                                                     const std::string& potential,
                                                     const std::string& overlap)
  : pkname_(pkname), tag_(tag), potential_(potential), overlap_(overlap)
{
  if (overlap_ == std::string("")) { overlap_ = potential_; }
};


void
UpwindPotentialDifference::Update(const CompositeVector& data,
                                  CompositeVector& uw_data,
                                  const State& S,
                                  const Teuchos::Ptr<Debugger>& db) const
{
  const CompositeVector& potential = S.Get<CompositeVector>(potential_, tag_);
  const CompositeVector& overlap = S.Get<CompositeVector>(overlap_, tag_);
  CalculateCoefficientsOnFaces(data, potential, overlap, uw_data);
};


void
UpwindPotentialDifference::CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                                        const CompositeVector& potential,
                                                        const CompositeVector& overlap,
                                                        CompositeVector& face_coef) const
{
  AMANZI_ASSERT(cell_coef.Ghosted());

  // initialize the cell coefficients
  if (face_coef.HasComponent("cell")) { face_coef.ViewComponent("cell", true)->PutScalar(1.0); }

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef.Mesh();
  std::vector<int> dirs;
  double eps = 1.e-16;

  // communicate ghosted cells
  cell_coef.ScatterMasterToGhosted("cell");
  potential.ScatterMasterToGhosted("cell");
  overlap.ScatterMasterToGhosted("cell");

  Epetra_MultiVector& face_coef_f = *face_coef.ViewComponent("face", false);
  const Epetra_MultiVector& overlap_c = *overlap.ViewComponent("cell", true);
  const Epetra_MultiVector& potential_c = *potential.ViewComponent("cell", true);
  Teuchos::RCP<const Epetra_MultiVector> potential_f;
  if (potential.HasComponent("face")) potential_f = potential.ViewComponent("face", false);
  const Epetra_MultiVector& cell_coef_c = *cell_coef.ViewComponent("cell", true);

  int nfaces = face_coef.size("face", false);
  for (unsigned int f = 0; f != nfaces; ++f) {
    auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);

    if (cells.size() == 1) {
      if (potential_f != Teuchos::null) {
        if (potential_c[0][cells[0]] >= (*potential_f)[0][f]) {
          face_coef_f[0][f] = cell_coef_c[0][cells[0]];
        }
      } else {
        face_coef_f[0][f] = cell_coef_c[0][cells[0]];
      }
    } else {
      // Determine the size of the overlap region, a smooth transition region
      // near zero potential difference.
      double ol0 = std::max(0., overlap_c[0][cells[0]]);
      double ol1 = std::max(0., overlap_c[0][cells[1]]);

      double flow_eps = 0.0;
      if ((ol0 > 0) || (ol1 > 0)) { flow_eps = (ol0 * ol1) / (ol0 + ol1); }
      flow_eps = std::max(flow_eps, eps);

      // Determine the coefficient.
      if (potential_c[0][cells[0]] - potential_c[0][cells[1]] > flow_eps) {
        face_coef_f[0][f] = cell_coef_c[0][cells[0]];
      } else if (potential_c[0][cells[1]] - potential_c[0][cells[0]] > flow_eps) {
        face_coef_f[0][f] = cell_coef_c[0][cells[1]];
      } else {
        // Parameterization of a linear scaling between upwind and downwind.
        double param;
        if (flow_eps < 2 * eps) {
          param = 0.5;
        } else {
          param = (potential_c[0][cells[1]] - potential_c[0][cells[0]]) / (2 * flow_eps) + 0.5;
        }
        AMANZI_ASSERT(param >= 0.0);
        AMANZI_ASSERT(param <= 1.0);
        face_coef_f[0][f] =
          cell_coef_c[0][cells[1]] * param + cell_coef_c[0][cells[0]] * (1. - param);
      }
    }
  }
};


void
UpwindPotentialDifference::UpdateDerivatives(
  const Teuchos::Ptr<State>& S,
  std::string potential_key,
  const CompositeVector& dconductivity,
  const std::vector<int>& bc_markers,
  const std::vector<double>& bc_values,
  std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>>* Jpp_faces) const
{
  double eps = 1.e-16;

  // Grab derivatives
  dconductivity.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& dcell_v = *dconductivity.ViewComponent("cell", true);

  AMANZI_ASSERT(dconductivity.Ghosted());

  // Grab potential
  AMANZI_ASSERT(potential_key == potential_);
  Teuchos::RCP<const CompositeVector> pres = S->GetPtr<CompositeVector>(potential_key, tag_);
  pres->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& pres_v = *pres->ViewComponent("cell", true);

  // Grab overlap
  Teuchos::RCP<const CompositeVector> overlap = S->GetPtr<CompositeVector>(overlap_, tag_);
  overlap->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& overlap_c = *overlap->ViewComponent("cell", true);

  // Grab mesh and allocate space
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = dconductivity.Mesh();
  unsigned int nfaces_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  Jpp_faces->resize(nfaces_owned);

  // workspace
  double dK_dp[2];
  double p[2];

  for (unsigned int f = 0; f != nfaces_owned; ++f) {
    auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int mcells = cells.size();

    // create the local matrix
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>> Jpp =
      Teuchos::rcp(new Teuchos::SerialDenseMatrix<int, double>(mcells, mcells));
    (*Jpp_faces)[f] = Jpp;

    if (mcells == 1) {
      if (bc_markers[f] == Operators::OPERATOR_BC_DIRICHLET) {
        // determine flux
        p[0] = pres_v[0][cells[0]];
        p[1] = bc_values[f];
        double dp = p[0] - p[1];

        if (p[0] > p[1]) {
          (*Jpp)(0, 0) = dp * mesh->getFaceArea(f) * dK_dp[0];
        } else {
          (*Jpp)(0, 0) = 0.;
        }
      } else {
        (*Jpp)(0, 0) = 0.;
      }

    } else {
      p[0] = pres_v[0][cells[0]];
      p[1] = pres_v[0][cells[1]];

      // Determine the size of the overlap region, a smooth transition region
      // near zero potential difference.
      double ol0 = std::max(0., overlap_c[0][cells[0]]);
      double ol1 = std::max(0., overlap_c[0][cells[1]]);

      double flow_eps = 0.0;
      if ((ol0 > 0) || (ol1 > 0)) { flow_eps = (ol0 * ol1) / (ol0 + ol1); }
      flow_eps = std::max(flow_eps, eps);

      // Determine the coefficient.
      if (p[0] - p[1] > flow_eps) {
        dK_dp[0] = dcell_v[0][cells[0]];
        dK_dp[1] = 0.;

      } else if (p[1] - p[0] > flow_eps) {
        dK_dp[0] = 0.;
        dK_dp[1] = dcell_v[0][cells[1]];

      } else {
        // Parameterization of a linear scaling between upwind and downwind.
        double param;
        if (flow_eps < 2 * eps) {
          param = 0.5;
        } else {
          param = (p[1] - p[0]) / (2 * flow_eps) + 0.5;
        }
        AMANZI_ASSERT(param >= 0.0);
        AMANZI_ASSERT(param <= 1.0);

        dK_dp[0] = (1. - param) * dcell_v[0][cells[0]];
        dK_dp[1] = param * dcell_v[0][cells[1]];
      }

      (*Jpp)(0, 0) = (p[0] - p[1]) * mesh->getFaceArea(f) * dK_dp[0];
      (*Jpp)(0, 1) = (p[0] - p[1]) * mesh->getFaceArea(f) * dK_dp[1];
      (*Jpp)(1, 0) = -(*Jpp)(0, 0);
      (*Jpp)(1, 1) = -(*Jpp)(0, 1);
    }
  }
}


} // namespace Operators
} // namespace Amanzi
