/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Op.hh"
#include "PDE_DiffusionWithGravity.hh"
#include "PDE_Accumulation.hh"

#include "pk_helpers.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (p + rho*g*z)
// -------------------------------------------------------------
void
Richards::ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  // force mass matrices to change
  if (!deform_key_.empty() &&
      S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " matrix"))
    matrix_diff_->SetTensorCoefficient(S_->GetPtr<TensorVector>(perm_key_, Tags::DEFAULT));

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(tag);

  // update the stiffness matrices
  matrix_->Zero();
  S_->GetEvaluator(mass_dens_key_, tag).Update(*S_, name_);
  matrix_diff_->SetDensity(S_->GetPtr<CompositeVector>(mass_dens_key_, tag));
  matrix_diff_->SetScalarCoefficient(S_->GetPtr<CompositeVector>(uw_coef_key_, tag), Teuchos::null);
  Teuchos::RCP<const CompositeVector> pres = S_->GetPtrW<CompositeVector>(key_, tag, name_);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());

  // derive fluxes
  Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(flux_key_, tag, name_);
  matrix_diff_->UpdateFlux(pres.ptr(), flux.ptr());
  PKHelpers::changedEvaluatorPrimary(flux_key_, tag, *S_);

  // apply BCs to the local matrices
  matrix_diff_->ApplyBCs(true, true, true);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres, *g);
};


// -------------------------------------------------------------
// Accumulation of water term du/dt
// -------------------------------------------------------------
void
Richards::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g)
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  // update the water content at both the old and new times.
  S_->GetEvaluator(conserved_key_, tag_next_).Update(*S_, name_);
  // S_->GetEvaluator(conserved_key_, tag_current_).Update(*S_, name_); // for the future...

  // get these fields
  Teuchos::RCP<const CompositeVector> wc1 = S_->GetPtr<CompositeVector>(conserved_key_, tag_next_);
  Teuchos::RCP<const CompositeVector> wc0 =
    S_->GetPtr<CompositeVector>(conserved_key_, tag_current_);

  // Water content only has cells, while the residual has cells and faces.
  g->getComponent("cell", false)
    ->update(1.0 / dt,
             *wc1->getComponent("cell", false),
             -1.0 / dt,
             *wc0->getComponent("cell", false),
             1.0);
};


// ---------------------------------------------------------------------
// Add in mass source, in units of mol / m^3 s
// ---------------------------------------------------------------------
void
Richards::AddSources_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // external sources of energy
  if (is_source_term_) {
    // Update the source term
    S_->GetEvaluator(source_key_, tag).Update(*S_, name_);
    g->getComponent("cell", false)
      ->elementWiseMultiply(
        -1.0,
        *S_->Get<CompositeVector>(source_key_, tag).getComponent("cell", false)->getVector(0),
        *S_->Get<CompositeVector>(cell_vol_key_, tag).getComponent("cell", false),
        1.0);

    db_->WriteVector("  source", S_->GetPtr<CompositeVector>(source_key_, tag).ptr(), false);
    db_->WriteVector("res (src)", g, false);
  }
}


void
Richards::AddSourcesToPrecon_(double h)
{
  // external sources of energy (temperature dependent source)
  if (is_source_term_ && !explicit_source_ && source_term_is_differentiable_ &&
      S_->GetEvaluator(source_key_, tag_next_).IsDifferentiableWRT(*S_, key_, tag_next_)) {
    S_->GetEvaluator(source_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
    preconditioner_acc_->AddAccumulationTerm(
      S_->GetDerivative<CompositeVector>(source_key_, tag_next_, key_, tag_next_),
      -1.0,
      "cell",
      true);
  }
}


// void
// Richards::UpdateVelocity_(const Tag& tag)
// {
//   AMANZI_ASSERT(tag == Tags::NEXT); // what else would this be?

//   const Epetra_MultiVector& flux =
//     *S_->Get<CompositeVector>(flux_key_, tag_next_).viewComponent("face", true);

//   S_->GetEvaluator(molar_dens_key_, tag_next_).Update(*S_, name_);
//   const Epetra_MultiVector& nliq_c =
//     *S_->Get<CompositeVector>(molar_dens_key_, tag_next_).viewComponent("cell", false);
//   Epetra_MultiVector& velocity =
//     *S_->GetW<CompositeVector>(velocity_key_, tag, name_).viewComponent("cell", true);

//   int d(mesh_->getSpaceDimension());
//   AmanziGeometry::Point local_velocity(d);

//   Teuchos::LAPACK<int, double> lapack;
//   Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
//   double rhs[d];

//   int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
//   AmanziMesh::Entity_ID_List faces;
//   for (int c = 0; c != ncells_owned; ++c) {
//     faces = mesh_->getCellFaces(c);
//     int nfaces = faces.size();

//     for (int i = 0; i != d; ++i) rhs[i] = 0.0;
//     matrix.putScalar(0.0);

//     for (int n = 0; n != nfaces; ++n) { // populate least-square matrix
//       int f = faces[n];
//       const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
//       double area = mesh_->getFaceArea(f);

//       for (int i = 0; i != d; ++i) {
//         rhs[i] += normal[i] * flux[0][f];
//         matrix(i, i) += normal[i] * normal[i];
//         for (int j = i + 1; j < d; ++j) { matrix(j, i) = matrix(i, j) += normal[i] * normal[j]; }
//       }
//     }

//     int info;
//     lapack.POSV('U', d, 1, matrix.values(), d, rhs, d, &info);

//     for (int i = 0; i != d; ++i) velocity[i][c] = rhs[i] / nliq_c[0][c];
//   }
// }


} // namespace Flow
} // namespace Amanzi
