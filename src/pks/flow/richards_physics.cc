/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Evaluator.hh"
#include "Op.hh"

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
    matrix_diff_->SetTensorCoefficient(K_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(tag);

  // update the matrix
  matrix_->Init();

  S_->GetEvaluator(mass_dens_key_, tag).Update(*S_, name_);
  matrix_diff_->SetDensity(S_->GetPtr<CompositeVector>(mass_dens_key_, tag));
  matrix_diff_->SetScalarCoefficient(S_->GetPtr<CompositeVector>(uw_coef_key_, tag), Teuchos::null);

  Teuchos::RCP<const CompositeVector> pres = S_->GetPtrW<CompositeVector>(key_, tag, name_);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // derive fluxes
  Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(flux_key_, tag, name_);
  matrix_diff_->UpdateFlux(pres.ptr(), flux.ptr());
  changedEvaluatorPrimary(flux_key_, tag, *S_);

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
  g->ViewComponent("cell", false)
    ->Update(1.0 / dt,
             *wc1->ViewComponent("cell", false),
             -1.0 / dt,
             *wc0->ViewComponent("cell", false),
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
    Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

    // Update the source term
    S_->GetEvaluator(source_key_, tag).Update(*S_, name_);
    const Epetra_MultiVector& source1 =
      *S_->Get<CompositeVector>(source_key_, tag).ViewComponent("cell", false);

    const Epetra_MultiVector& cv =
      *S_->Get<CompositeVector>(Keys::getKey(domain_, "cell_volume"), tag)
         .ViewComponent("cell", false);

    // Add into residual
    unsigned int ncells = g_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) { g_c[0][c] -= source1[0][c] * cv[0][c]; }

    db_->WriteVector("  source", S_->GetPtr<CompositeVector>(source_key_, tag).ptr(), false);
    db_->WriteVector("res (src)", g, false);
  }
}


void
Richards::AddSourcesToPrecon_(double h)
{
  // nonlinear water sources include a dQ/dp term into PC
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


// -------------------------------------------------------------
// Convert abs perm vector to tensor.
// -------------------------------------------------------------
void
Richards::SetAbsolutePermeabilityTensor_(const Tag& tag)
{
  // currently assumes isotropic perm, should be updated
  S_->GetEvaluator(perm_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& perm =
    *S_->Get<CompositeVector>(perm_key_, tag).ViewComponent("cell", false);
  unsigned int ncells = perm.MyLength();
  unsigned int ndofs = perm.NumVectors();
  int space_dim = mesh_->space_dimension();
  if (ndofs == 1) { // isotropic
    for (unsigned int c = 0; c != ncells; ++c) { (*K_)[c](0, 0) = perm[0][c] * perm_scale_; }
  } else if (ndofs == 2 && space_dim == 3) {
    // horizontal and vertical perms
    for (unsigned int c = 0; c != ncells; ++c) {
      (*K_)[c](0, 0) = perm[0][c] * perm_scale_;
      (*K_)[c](1, 1) = perm[0][c] * perm_scale_;
      (*K_)[c](2, 2) = perm[1][c] * perm_scale_;
    }
  } else if (ndofs >= space_dim) {
    // diagonal tensor
    for (unsigned int dim = 0; dim != space_dim; ++dim) {
      for (unsigned int c = 0; c != ncells; ++c) {
        (*K_)[c](dim, dim) = perm[dim][c] * perm_scale_;
      }
    }
    if (ndofs > space_dim) {
      // full tensor
      if (ndofs == 3) { // 2D
        for (unsigned int c = 0; c != ncells; ++c) {
          (*K_)[c](0, 1) = (*K_)[c](1, 0) = perm[2][c] * perm_scale_;
        }
      } else if (ndofs == 6) { // 3D
        for (unsigned int c = 0; c != ncells; ++c) {
          (*K_)[c](0, 1) = (*K_)[c](1, 0) = perm[3][c] * perm_scale_; // xy & yx
          (*K_)[c](0, 2) = (*K_)[c](2, 0) = perm[4][c] * perm_scale_; // xz & zx
          (*K_)[c](1, 2) = (*K_)[c](2, 1) = perm[5][c] * perm_scale_; // yz & zy
        }
      }
    }
  } else {
    // ERROR -- unknown perm type
    AMANZI_ASSERT(0);
  }
};


void
Richards::UpdateVelocity_(const Tag& tag)
{
  AMANZI_ASSERT(tag == Tags::NEXT); // what else would this be?

  const Epetra_MultiVector& flux =
    *S_->Get<CompositeVector>(flux_key_, tag_next_).ViewComponent("face", true);

  S_->GetEvaluator(molar_dens_key_, tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& nliq_c =
    *S_->Get<CompositeVector>(molar_dens_key_, tag_next_).ViewComponent("cell", false);
  Epetra_MultiVector& velocity =
    *S_->GetW<CompositeVector>(velocity_key_, tag, name_).ViewComponent("cell", true);

  int d(mesh_->space_dimension());
  AmanziGeometry::Point local_velocity(d);

  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
  double rhs[d];

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int i = 0; i != d; ++i) rhs[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n = 0; n != nfaces; ++n) { // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i = 0; i != d; ++i) {
        rhs[i] += normal[i] * flux[0][f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i + 1; j < d; ++j) { matrix(j, i) = matrix(i, j) += normal[i] * normal[j]; }
      }
    }

    int info;
    lapack.POSV('U', d, 1, matrix.values(), d, rhs, d, &info);

    for (int i = 0; i != d; ++i) velocity[i][c] = rhs[i] / nliq_c[0][c];
  }
}


// // -------------------------------------------------------------
// // Diffusion term, div -\phi s \tau n D grad \omega
// // -------------------------------------------------------------
// void Richards::AddVaporDiffusionResidual_(const Teuchos::Ptr<State>& S,
//         const Teuchos::Ptr<CompositeVector>& g) {

//   //res_vapor = Teuchos::rcp(new CompositeVector(*S_->GetPtr<CompositeVector>("pressure")));
//   //res_vapor = Teuchos::rcp(new CompositeVector(*g));
//   res_vapor->PutScalar(0.0);
//   //Teuchos::RCP<CompositeVector> res_en = Teuchos::rcp(new CompositeVector(*g));
//   //res_en->PutScalar(0.);

//   // derive fluxes
//   Teuchos::RCP<const CompositeVector> pres   = S_->GetPtr<CompositeVector>("pressure");
//   Teuchos::RCP<const CompositeVector> temp   = S_->GetPtr<CompositeVector>("temperature");


//   Teuchos::RCP<CompositeVector> vapor_diff_pres = S_->GetPtrW<CompositeVector>("vapor_diffusion_pressure", name_);
//   Teuchos::RCP<CompositeVector> vapor_diff_temp = S_->GetPtrW<CompositeVector>("vapor_diffusion_temperature", name_);

//   ///****** Compute contribution for pressure gradient

//   //Epetra_MultiVector& coef_pr = *vapor_diff_pres->ViewComponent("cell",false);
//   ComputeVaporDiffusionCoef(S, vapor_diff_pres, "pressure");

//   // update the stiffness matrix
//   matrix_vapor_->CreateMFDstiffnessMatrices(vapor_diff_pres.ptr());
//   matrix_vapor_->CreateMFDrhsVectors();
//   // assemble the stiffness matrix
//   //matrix_vapor_->ApplyBoundaryConditions(bc_markers_, bc_values_, false);
//   //  matrix_vapor_->AssembleGlobalMatrices();
//   // calculate the residual
//   matrix_vapor_->ComputeNegativeResidual(*pres, res_vapor.ptr());

//   g->Update(1., *res_vapor, 1.);

//   res_vapor->PutScalar(0.0);

//   ///****** Compute contribution for temperature gradient
//   Epetra_MultiVector& coef_tm = *S_->GetPtrW<CompositeVector>("vapor_diffusion_temperature", name_)
//                                    ->ViewComponent("cell",false);

//   ComputeVaporDiffusionCoef(S, vapor_diff_temp, "temperature");
//   // update the stiffness matrix
//   matrix_vapor_->CreateMFDstiffnessMatrices(vapor_diff_temp.ptr());
//   matrix_vapor_->CreateMFDrhsVectors();
//   // assemble the stiffness matrix
//   //matrix_vapor_->ApplyBoundaryConditions(bc_markers_, bc_values_, false);
//   //  matrix_vapor_->AssembleGlobalMatrices();
//   // calculate the residual
//   matrix_vapor_->ComputeNegativeResidual(*temp, res_vapor.ptr());

//   g->Update(1., *res_vapor, 1.);


// }

//   void Richards::ComputeVaporDiffusionCoef(const Teuchos::Ptr<State>& S,
//                                           Teuchos::RCP<CompositeVector>& vapor_diff,
//                                           std::string var_name){

//    Epetra_MultiVector& diff_coef = *vapor_diff->ViewComponent("cell",false);

//    S_->GetEvaluator("molar_density_liquid").Update(S.ptr(), name_);
//    const Epetra_MultiVector& n_l = *S_->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell",false);

//    S_->GetEvaluator("molar_density_gas").Update(S.ptr(), name_);
//    const Epetra_MultiVector& n_g = *S_->Get<CompositeVector>("molar_density_gas").ViewComponent("cell",false);

//    S_->GetEvaluator("porosity").Update(S.ptr(), name_);
//    const Epetra_MultiVector& phi = *S_->Get<CompositeVector>("porosity").ViewComponent("cell",false);

//    S_->GetEvaluator("saturation_gas").Update(S.ptr(), name_);
//    const Epetra_MultiVector& s_g = *S_->Get<CompositeVector>("saturation_gas").ViewComponent("cell",false);

//    S_->GetEvaluator("mol_frac_gas").Update(S.ptr(), name_);
//    const Epetra_MultiVector& mlf_g = *S_->Get<CompositeVector>("mol_frac_gas").ViewComponent("cell",false);

//    std::string key_t = "temperature";
//    S_->GetEvaluator("mol_frac_gas")->HasFieldDerivativeChanged(S.ptr(), name_, key_t);
//    const Epetra_MultiVector& dmlf_g_dt = *S_->Get<CompositeVector>("dmol_frac_gas_dtemperature").ViewComponent("cell",false);

//    const Epetra_MultiVector& temp = *S_->Get<CompositeVector>("temperature").ViewComponent("cell",false);
//    const Epetra_MultiVector& pressure = *S_->Get<CompositeVector>("pressure").ViewComponent("cell",false);
//    const double& Patm = *S_->GetScalarData("atmospheric_pressure", Tags::DEFAULT);
//    const double R = 8.3144621;

//    unsigned int ncells = diff_coef.MyLength();

//    const double a = 4./3.;
//    const double b = 10./3.;
//    const double D_ref = 0.282;
//    const double P_ref = Patm;
//    const double T_ref = 298;
//    double D;

//    for (unsigned int c=0; c!=ncells; ++c){

//      D = D_ref*(P_ref/Patm)*pow(temp[0][c]/T_ref, 1.8);

//      diff_coef[0][c] = D*pow(phi[0][c], a)*pow(s_g[0][c], b)*n_g[0][c];
//      diff_coef[0][c] *= exp(-(Patm - pressure[0][c])/(n_l[0][c]*R*temp[0][c]));
//    }

//    if (var_name == "pressure"){
//      //cout<<"Pressure vapor_diff\n";
//      for (unsigned int c=0; c!=ncells; ++c){
//        diff_coef[0][c] *= mlf_g[0][c] * (1./ (n_l[0][c]*R*temp[0][c]));
//        //diff_coef[0][c] *= 0.;
//        //cout<<diff_coef[0][c]<<" ";
//      }
//      //cout<<endl;
//    }
//    else if (var_name == "temperature"){
//      for (unsigned int c=0; c!=ncells; ++c){
//        diff_coef[0][c] *= (1./Patm)*dmlf_g_dt[0][c] + mlf_g[0][c]* (Patm - pressure[0][c])/ (n_l[0][c]*R*temp[0][c]*temp[0][c]);
//        //diff_coef[0][c] =0;
//      }

//    }
//    else{
//      // Unknown variable name
//      AMANZI_ASSERT(0);
//    }

// }

} // namespace Flow
} // namespace Amanzi
