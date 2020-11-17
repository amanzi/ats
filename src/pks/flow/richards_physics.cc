/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Neil Carlson (version 1)
  Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "FieldEvaluator.hh"
#include "Op.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad (p + rho*g*z)
// -------------------------------------------------------------
void Richards::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  // update the rel perm according to the scheme of choice
  bool update = UpdatePermeabilityData_(S.ptr());

  // update the matrix
  matrix_->Init();

  S->GetFieldEvaluator(mass_dens_key_)->HasFieldChanged(S, name_);
  matrix_diff_->SetDensity(S->GetFieldData(mass_dens_key_));
  matrix_diff_->SetScalarCoefficient(S->GetFieldData(uw_coef_key_), Teuchos::null);

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_, name_);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // derive fluxes
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData(flux_key_, name_);
  matrix_diff_->UpdateFlux(pres.ptr(), flux.ptr());
  if (S == S_next_.ptr()) flux_pvfe_->SetFieldAsChanged(S);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres, *g);
};


// -------------------------------------------------------------
// Accumulation of water term du/dt
// -------------------------------------------------------------
void Richards::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // update the water content at both the old and new times.
  S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_inter_.ptr(), name_);
    
  // get these fields
  Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData(conserved_key_);
  Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData(conserved_key_);

  // Water content only has cells, while the residual has cells and faces.
  g->ViewComponent("cell",false)->Update(1.0/dt, *wc1->ViewComponent("cell",false),
          -1.0/dt, *wc0->ViewComponent("cell",false), 1.0);
  
  db_->WriteVector("res (acc)", g, true);

};


// ---------------------------------------------------------------------
// Add in mass source, in units of mol / m^3 s
// ---------------------------------------------------------------------
void Richards::AddSources_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // external sources of energy
  if (is_source_term_) {
    Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

    // Update the source term
    S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
    const Epetra_MultiVector& source1 =
        *S->GetFieldData(source_key_)->ViewComponent("cell",false);

    const Epetra_MultiVector& cv =
      *S->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

    // Add into residual
    unsigned int ncells = g_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      g_c[0][c] -= source1[0][c] * cv[0][c];
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Adding external source term" << std::endl;
      db_->WriteVector("  Q_ext", S->GetFieldData(source_key_).ptr(), false);
    }  
    db_->WriteVector("res (src)", g, false);
  }
}


void Richards::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
  // external sources of energy (temperature dependent source)
  if (is_source_term_ && !explicit_source_ && source_term_is_differentiable_ &&
      S->GetFieldEvaluator(source_key_)->IsDependency(S, key_)) {

    S->GetFieldEvaluator(source_key_)->HasFieldDerivativeChanged(S, name_, key_);
    Key dsource_dp_key = Keys::getDerivKey(source_key_, key_);

    preconditioner_acc_->AddAccumulationTerm(*S->GetFieldData(dsource_dp_key), -1.0, "cell", true);
  }
}


// -------------------------------------------------------------
// Convert abs perm vector to tensor.
// -------------------------------------------------------------
void Richards::SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S) {
  // currently assumes isotropic perm, should be updated
  S->GetFieldEvaluator(perm_key_)->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& perm = *S->GetFieldData(perm_key_)
      ->ViewComponent("cell",false);
  unsigned int ncells = perm.MyLength();
  unsigned int ndofs = perm.NumVectors();

  if (ndofs == 1) { // isotropic
    for (unsigned int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = perm[0][c] * perm_scale_;
    }
  } else if (ndofs == 2 && S->GetMesh()->space_dimension() == 3) {
    // horizontal and vertical perms
    for (unsigned int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = perm[0][c] * perm_scale_;
      (*K_)[c](1, 1) = perm[0][c] * perm_scale_;
      (*K_)[c](2, 2) = perm[1][c] * perm_scale_;
    }
  } else if (ndofs >= S->GetMesh()->space_dimension()) {
    // diagonal tensor
    for (unsigned int dim=0; dim!=S->GetMesh()->space_dimension(); ++dim) {
      for (unsigned int c=0; c!=ncells; ++c) {
        (*K_)[c](dim, dim) = perm[dim][c] * perm_scale_;
      }
    }
      if (ndofs > S->GetMesh()->space_dimension()) {
      // full tensor
      if (ndofs == 3) { // 2D
        for (unsigned int c=0; c!=ncells; ++c) {
          (*K_)[c](0,1) = (*K_)[c](1,0) = perm[2][c] * perm_scale_;
        }
      } else if (ndofs == 6) { // 3D
        for (unsigned int c=0; c!=ncells; ++c) {
          (*K_)[c](0,1) = (*K_)[c](1,0) = perm[3][c] * perm_scale_; // xy & yx
          (*K_)[c](0,2) = (*K_)[c](2,0) = perm[4][c] * perm_scale_; // xz & zx
          (*K_)[c](1,2) = (*K_)[c](2,1) = perm[5][c] * perm_scale_; // yz & zy
        }
      }
    }
  } else {
    // ERROR -- unknown perm type
    AMANZI_ASSERT(0);
  }
};


void
Richards::UpdateVelocity_(const Teuchos::Ptr<State>& S) {
  const Epetra_MultiVector& flux = *S->GetFieldData(flux_key_)
      ->ViewComponent("face", true);

  S->GetFieldEvaluator(molar_dens_key_)->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& nliq_c = *S->GetFieldData(molar_dens_key_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& velocity = *S->GetFieldData(velocity_key_, name_)
      ->ViewComponent("cell", true);

  int d(mesh_->space_dimension());
  AmanziGeometry::Point local_velocity(d);

  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
  double rhs[d];

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int i=0; i!=d; ++i) rhs[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n=0; n!=nfaces; ++n) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i=0; i!=d; ++i) {
        rhs[i] += normal[i] * flux[0][f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i+1; j < d; ++j) {
          matrix(j, i) = matrix(i, j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', d, 1, matrix.values(), d, rhs, d, &info);

    for (int i=0; i!=d; ++i) velocity[i][c] = rhs[i] / nliq_c[0][c];
  }
}


// // -------------------------------------------------------------
// // Diffusion term, div -\phi s \tau n D grad \omega
// // -------------------------------------------------------------
// void Richards::AddVaporDiffusionResidual_(const Teuchos::Ptr<State>& S,
//         const Teuchos::Ptr<CompositeVector>& g) {

//   //res_vapor = Teuchos::rcp(new CompositeVector(*S->GetFieldData("pressure"))); 
//   //res_vapor = Teuchos::rcp(new CompositeVector(*g)); 
//   res_vapor->PutScalar(0.0);
//   //Teuchos::RCP<CompositeVector> res_en = Teuchos::rcp(new CompositeVector(*g)); 
//   //res_en->PutScalar(0.);

//   // derive fluxes
//   Teuchos::RCP<const CompositeVector> pres   = S->GetFieldData("pressure");
//   Teuchos::RCP<const CompositeVector> temp   = S->GetFieldData("temperature");


//   Teuchos::RCP<CompositeVector> vapor_diff_pres = S->GetFieldData("vapor_diffusion_pressure", name_);
//   Teuchos::RCP<CompositeVector> vapor_diff_temp = S->GetFieldData("vapor_diffusion_temperature", name_);

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
//   Epetra_MultiVector& coef_tm = *S->GetFieldData("vapor_diffusion_temperature", name_)
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

//    S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S.ptr(), name_);
//    const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);

//    S->GetFieldEvaluator("molar_density_gas")->HasFieldChanged(S.ptr(), name_);
//    const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);

//    S->GetFieldEvaluator("porosity")->HasFieldChanged(S.ptr(), name_);
//    const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);

//    S->GetFieldEvaluator("saturation_gas")->HasFieldChanged(S.ptr(), name_);
//    const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);

//    S->GetFieldEvaluator("mol_frac_gas")->HasFieldChanged(S.ptr(), name_);
//    const Epetra_MultiVector& mlf_g = *S->GetFieldData("mol_frac_gas")->ViewComponent("cell",false);

//    std::string key_t = "temperature";
//    S->GetFieldEvaluator("mol_frac_gas")->HasFieldDerivativeChanged(S.ptr(), name_, key_t);
//    const Epetra_MultiVector& dmlf_g_dt = *S->GetFieldData("dmol_frac_gas_dtemperature")->ViewComponent("cell",false);

//    const Epetra_MultiVector& temp = *S->GetFieldData("temperature")->ViewComponent("cell",false);
//    const Epetra_MultiVector& pressure = *S->GetFieldData("pressure")->ViewComponent("cell",false);
//    const double& Patm = *S->GetScalarData("atmospheric_pressure");
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

} //namespace
} //namespace
