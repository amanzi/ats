/* -*-  mode++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Svetlana Tokareva

Solves:

c rho dT/dt + q dot grad h = div Ke grad T + S
------------------------------------------------------------------------- */

#include "advection.hh"
#include "FieldEvaluator.hh"
#include "lake_thermo_pk.hh"
#include "Op.hh"
#include "pk_helpers.hh"

namespace Amanzi {
namespace LakeThermo {

// -------------------------------------------------------------
// Accumulation of energy term c rho dT/dt
// -------------------------------------------------------------
void Lake_Thermo_PK::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

//  // update the temperature at both the old and new times.
//  S_next_->GetFieldEvaluator(temperature_key_)->HasFieldChanged(S_next_.ptr(), name_);
//  S_inter_->GetFieldEvaluator(temperature_key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // get the temperature at each time
  Teuchos::RCP<const CompositeVector> T1 = S_next_->GetFieldData(temperature_key_);
  Teuchos::RCP<const CompositeVector> T0 = S_inter_->GetFieldData(temperature_key_);

  // evaluate density
  S_inter_->GetFieldEvaluator(density_key_)->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& rho =
      *S_inter_->GetFieldData(density_key_)->ViewComponent("cell",false);

  // evaluate heat capacity
  S_inter_->GetFieldEvaluator(heat_capacity_key_)->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& cp =
      *S_inter_->GetFieldData(heat_capacity_key_)->ViewComponent("cell",false);

  // Update the residual with the accumulation of energy over the

  const Epetra_MultiVector& T1_c = *S_next_->GetFieldData(temperature_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& T0_c = *S_inter_->GetFieldData(temperature_key_)->ViewComponent("cell", false);

  const Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

  const Epetra_MultiVector& cv =
        *S_inter_->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

  // get the energy at each time

  S_next_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_inter_.ptr(), name_);  

  Teuchos::RCP<const CompositeVector> e1 = S_next_->GetFieldData(energy_key_);
  Teuchos::RCP<const CompositeVector> e0 = S_inter_->GetFieldData(energy_key_);   

  const Epetra_MultiVector& e1_c = *S_next_->GetFieldData(energy_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& e0_c = *S_inter_->GetFieldData(energy_key_)->ViewComponent("cell", false); 

  unsigned int ncells = g_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    g_c[0][c] += (e1_c[0][c]-e0_c[0][c])/dt * cv[0][c];
  }    

  // unsigned int ncells = g_c.MyLength();
  // for (unsigned int c=0; c!=ncells; ++c) {
  //   g_c[0][c] += cp[0][c]*rho[0][c]*(T1_c[0][c]-T0_c[0][c])/dt * cv[0][c] ;
  // }

};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, q dot grad h.
// -------------------------------------------------------------
void Lake_Thermo_PK::AddAdvection_(const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& g, bool negate) {

  // set up the operator

  // NOTE: fluxes are a MOLAR flux by choice of the flow pk, i.e.
  // [flux] =  mol/s

  // NOTE: this will be the eventual way to ensure it is up to date,
  // but there is no FieldEvaluator for darcy flux yet.  When there
  // is, we can take the evaluation out of Flow::commit_state(),
  // but for now we'll leave it there and assume it has been updated. --etc
  //  S->GetFieldEvaluator(flux_key_)->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_key_);

  // evaluate density
  S->GetFieldEvaluator(density_key_)->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& rho_v =
      *S->GetFieldData(density_key_)->ViewComponent("cell",false);

  // evaluate heat capacity
  S->GetFieldEvaluator(heat_capacity_key_)->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& cp_v =
      *S->GetFieldData(heat_capacity_key_)->ViewComponent("cell",false);

  Teuchos::ParameterList& param_list = plist_->sublist("met data");
  FunctionFactory fac;
  Teuchos::RCP<Function> r_func_ = Teuchos::rcp(fac.Create(param_list.sublist("precipitation")));
  Teuchos::RCP<Function> E_func_ = Teuchos::rcp(fac.Create(param_list.sublist("evaporation")));

  std::vector<double> args(1);
  args[0] = S->time();
  r_ = (*r_func_)(args);
  E_ = (*E_func_)(args);

  // Compute evaporartion rate
  S->GetFieldEvaluator(evaporation_rate_key_)->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& E_v =
      *S->GetFieldData(evaporation_rate_key_)->ViewComponent("cell",false);
  E_ = E_v[0][0]; // same everywhere

  // set precipitation and evaporation to zero in the winter (for now because we don't have a snow layer anyway)
  unsigned int ncells = g->ViewComponent("cell",false)->MyLength();
  // get temperature
  const Epetra_MultiVector& temp1 = *S->GetFieldData(temperature_key_)
              ->ViewComponent("cell",false);

  // set precipitation and evaporation to zero in the winter (for now because we don't have a snow layer anyway)
  // r_ = (temp1[0][ncells-1] < 273.15) ? 0. : r_;
  // E_ = (temp1[0][ncells-1] < 273.15) ? 0. : E_;

  std::ofstream outfile;

  // outfile.open("Precipitation.txt", std::ios_base::app); // append instead of overwrite
  // outfile.precision(5);
  // outfile.setf(std::ios::scientific,std::ios::floatfield);
  // outfile << S->time() << " " << r_ << "\n";

  double dhdt = r_ - E_ - R_s_ - R_b_;
  double B_w  = r_ - E_;
  
  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face", false);

  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  double dt = S_next_->time() - S_inter_->time();

  const Epetra_MultiVector& cv =
          *S_inter_->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

//    AmanziGeometry::Point normal = mesh_->face_normal(f);
//    normal /= norm(normal);

    double cp;     // cell-based but we need face values for the flux
    double rho;

    AmanziMesh::Entity_ID_List f_cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &f_cells);

    int orientation;
    AmanziGeometry::Point normal = mesh_->face_normal(f, false, f_cells[0], &orientation);

  //  double *area;
  //  AmanziGeometry::Point centroid;
  //  std::vector<AmanziGeometry::Point> normals;
  //  mesh_->compute_face_geometry_(f,area,centroid,normals);

   AmanziMesh::Entity_ID_List cellfaceids;
   std::vector<int> cellfacedirs;
   int dir = 1;

   if (f_cells.size() == 1) mesh_->cell_get_faces_and_dirs(f_cells[0], &cellfaceids, &cellfacedirs);
   if (f_cells.size() == 2) mesh_->cell_get_faces_and_dirs(f_cells[1], &cellfaceids, &cellfacedirs);

   for (int j = 0; j < cellfaceids.size(); j++) {
     if (cellfaceids[j] == f) {
       dir = cellfacedirs[j];
       break;
     }
   }

  //  normal = normal*dir;

    if (f_cells.size() == 1) {
      // boundary face, use the cell value
      cp  = cp_v[0][f_cells[0]];
      rho = rho_v[0][f_cells[0]];
    } else {
      AMANZI_ASSERT(f_cells.size() == 2);
      // interpolate between cells
      cp  = (cp_v[0][f_cells[0]]  + cp_v[0][f_cells[1]])  / 2.;
      rho = (rho_v[0][f_cells[0]] + rho_v[0][f_cells[1]]) / 2.;
    }

    double dhdt_c = dhdt;
    double B_w_c = B_w;

    dhdt_c = (temp1[0][ncells-1] < 273.15) ? 0. : dhdt;
    B_w_c = (temp1[0][ncells-1] < 273.15) ? 0. : B_w;

    flux_f[0][f] = -1.*cp*rho*(dhdt_c*(1.-xcf[2]) - B_w_c)/(h_+1.e-12); // / cv[0][f_cells[0]]; //* normal[2]; // / cv[0][f_cells[0]]; // * normal[2] ; // / cv[0][f_cells[0]]; //*normal[2]; // *normal or not?
    // flux_f[0][f] = -1.*(dhdt_c*(1.-xcf[2]) - B_w_c)/(h_+1.e-12) / cv[0][f_cells[0]];
  //  std::cout << "f = " << f << ", normal = " << normal << ", dir = " << dir << std::endl;
  //  for (int j = 0; j < f_cells.size(); j++) std::cout << "f_cells[" << j <<"] = " << f_cells[j] << std::endl;
  }
//  exit(0);

  db_->WriteVector(" adv flux", flux.ptr(), true);
  matrix_adv_->global_operator()->Init();
  matrix_adv_->Setup(*flux);
  matrix_adv_->SetBCs(bc_adv_, bc_adv_);
  matrix_adv_->UpdateMatrices(flux.ptr());

  // apply to temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);
  ApplyDirichletBCsToTemperature_(S.ptr());
  // ApplyDirichletBCsToEnergy_(S.ptr());
  matrix_adv_->ApplyBCs(false, true, false);

  // apply
  matrix_adv_->global_operator()->ComputeNegativeResidual(*temp, *g, false);

  // S->GetFieldEvaluator(energy_key_)->HasFieldChanged(S.ptr(), name_);
  // Teuchos::RCP<const CompositeVector> enrg = S->GetFieldData(energy_key_);
  // matrix_adv_->global_operator()->ComputeNegativeResidual(*enrg, *g, false); 
}

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void Lake_Thermo_PK::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& g) {
  // update the thermal conductivity
  UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  const Epetra_MultiVector& uw_cond_c = *S_next_->GetFieldData(uw_conductivity_key_)->ViewComponent("face", false);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  const Epetra_MultiVector& cv =
      *S->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

  // get temperature
  const Epetra_MultiVector& temp_v = *S->GetFieldData(temperature_key_)
            ->ViewComponent("cell",false);

  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);    

  for (int f = 0; f < nfaces_owned; f++) {

    uw_cond_c[0][f] = uw_cond_c[0][f]; ///cv[0][0]; ///cv[0][0]; // ///(h_*h_) ;///cv[0][0];

  }

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);

  // update the stiffness matrix
  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // calculate the residual
  matrix_diff_->global_operator()->ComputeNegativeResidual(*temp, *g);

};


// ---------------------------------------------------------------------
// Add in energy source, which are accumulated by a single evaluator.
// ---------------------------------------------------------------------
void Lake_Thermo_PK::AddSources_(const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  is_source_term_ = true;

  // external sources of energy
  if (is_source_term_) {
    Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

    // Update the source term
    //    S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
    //    const Epetra_MultiVector& source1 =
    //        *S->GetFieldData(source_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& cv =
        *S->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

    Teuchos::ParameterList& param_list = plist_->sublist("met data");
    FunctionFactory fac;
    Teuchos::RCP<Function> r_func_  = Teuchos::rcp(fac.Create(param_list.sublist("precipitation")));
    Teuchos::RCP<Function> E_func_  = Teuchos::rcp(fac.Create(param_list.sublist("evaporation")));
    Teuchos::RCP<Function> SS_func_ = Teuchos::rcp(fac.Create(param_list.sublist("solar radiation")));

    std::vector<double> args(1);
    args[0] = S->time();
    r_ = (*r_func_)(args);
    E_ = (*E_func_)(args);

    double SS = (*SS_func_)(args);

    // Compute evaporartion rate
    S->GetFieldEvaluator(evaporation_rate_key_)->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& E_v =
        *S->GetFieldData(evaporation_rate_key_)->ViewComponent("cell",false);
    E_ = E_v[0][0]; // same everywhere

    double dhdt = r_ - E_ - R_s_ - R_b_;
    double B_w  = r_ - E_;

    S->GetFieldEvaluator(density_key_)->HasFieldChanged(S.ptr(), name_);

    // evaluate density
    const Epetra_MultiVector& rho =
        *S->GetFieldData(density_key_)->ViewComponent("cell",false);

    // get temperature
    const Epetra_MultiVector& temp = *S->GetFieldData(temperature_key_)
              ->ViewComponent("cell",false);

    // get energy
    S->GetFieldEvaluator(energy_key_)->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& enrg = *S->GetFieldData(energy_key_)
              ->ViewComponent("cell",false);          

    // get conductivity
    S->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& lambda_c =
        *S->GetFieldData(conductivity_key_)->ViewComponent("cell",false);

    // evaluate heat capacity
    S->GetFieldEvaluator(heat_capacity_key_)->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& cp =
        *S->GetFieldData(heat_capacity_key_)->ViewComponent("cell",false);

    // water extinction coefficient
    alpha_e_w_ = 1.04;
    // ice extinction coefficient
    alpha_e_i_ = 1.0;

    double albedo_w = 0.06; // water albedo
    double albedo_i = 0.40; // ice albedo

    double alpha;

    // Add into residual
    unsigned int ncells = g_c.MyLength();

    // // set precipitation and evaporation to zero in the winter (for now because we don't have a snow layer anyway)
    // r_ = (temp[0][ncells-1] < 273.15) ? 0. : r_;
    // E_ = (temp[0][ncells-1] < 273.15) ? 0. : E_;

    for (unsigned int c=0; c!=ncells; ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

      alpha = (temp[0][c] < 273.15) ? alpha_e_i_ : alpha_e_w_;

      double hour_sec = 60.*60;
      double interval = 24.;
    //  SS = SS/(hour_sec*interval);    

     int n = int(S->time());
     int day = n / (24 * 3600);
     n = n % (24 * 3600);
     int hour = n / 3600;
     n %= 3600;
     int minutes = n / 60 ;
     n %= 60;
     int seconds = n;

     double pi = 3.1415;
     // FoxDen, Atqasuk
     double longitude = -164.456696;
     double latitude = 66.558765;
     double phi   = latitude*pi/180.0;
    //  double delta = longitude*pi/180.0; 
    double delta = 23.5*pi/180.0*cos(2*pi*(double(day)-173.0)/365.0);
     double theta = pi*(hour-12.0)/12.0;
     double sinh0 = sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(theta);
     sinh0 = fmax(sinh0,0.0); 

    //  SS *= sqrt(1.-sinh0*sinh0);

      // S0_ = SS;
      // S0_ = 0.95*SS; // FoxDen
      S0_ = 0.9*SS; // Atqasuk

      double albedo = (temp[0][ncells-1] < 273.15) ? albedo_i : albedo_w;

      S0_ *= (1.-albedo);

      double S0_w_ice = S0_*exp(-alpha_e_i_*h_ice_);
      double dhdt_c = dhdt;
      double B_w_c = B_w;  

      S0_ = (temp[0][ncells-1] < 273.15) ? S0_w_ice : S0_;
      dhdt_c = (temp[0][ncells-1] < 273.15) ? 0. : dhdt;
      B_w_c = (temp[0][ncells-1] < 273.15) ? 0. : B_w;

      double snow_depth = 0.;
      if (coupled_to_snow_) {
        const Epetra_MultiVector& snow_depth_v = *S->GetFieldData("snow-depth")->ViewComponent("cell",false);
        snow_depth = snow_depth_v[0][0] ;
      }
      if (snow_depth > 2.e-2) {
        S0_ = 0.;
      }

      // S0_ = (temp[0][ncells-1] < 273.15) ? 0. : S0_;

      // -1.* because I switched to vertical xi coordinate
      g_c[0][c] += -1.*( (S0_*exp(-alpha*h_*(1.-xc[2]))*(alpha*h_)/(h_+1.e-12) - cp[0][c]*rho[0][c]*temp[0][c]*dhdt_c/(h_+1.e-12)) ) * cv[0][c] ; // * cv[0][c] ;
      // g_c[0][c] += -1.*( (S0_*exp(-alpha*h_*(1.-xc[2]))*(alpha*h_)/(h_+1.e-12) - enrg[0][c]*dhdt_c/(h_+1.e-12)) ); // * cv[0][c] ;
      // g_c[0][c] += -1.*( (S0_*exp(-alpha*h_*(1.-xc[2]))*(alpha*h_)/(h_+1.e-12)) ) * cv[0][c] ;


      // g_c[0][c] += -1.*( (S0_*exp(-alpha*h_*(1.-xc[2]))*(alpha*h_)/(h_+1.e-6) - enrg[0][c]*dhdt_c/(h_+1.e-6)) ) * cv[0][c] ;

      /* TESTING
//      // manufactured solution 1: linear temperature distribution
//      double T0 = 280., T1 = 400.;
//      g_c[0][c] += cp_*rho[0][c]*(T1-T0)/(h_+1.e-6)*(dhdt*xc[2] - B_w) * cv[0][c]; // - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
      // manufactured solution 2: exponential temperature distribution
      double T0 = 280., T1 = 400.;
      double coef = log(T1/T0);
//      std::cout << "lambda_c[0][c] = " << lambda_c[0][c] << std::endl;
//      std::cout << "coef = " << coef << std::endl;
//      std::cout << "exp(coef*xc[2]) = " << exp(coef*xc[2]) << std::endl;
//      std::cout << "h_ = " << h_ << std::endl;
      g_c[0][c] += T0*coef*exp(coef*xc[2])*( 1./(h_+1.e-6)*lambda_c[0][c]*coef - cp_*rho[0][c]*(dhdt*xc[2] - B_w) )/(h_+1.e-6) * cv[0][c];
//          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      // manufactured solution 3: quadratic temperature distribution
//      double T0 = 280., T1 = 400.;
//      double zm = 0.5, Tm = 300.;
//      double d = T0;
//      double b = (Tm-T0-zm*zm*(T1-T0))/(zm*(1.-zm));
//      double a = T1-T0-b;
////      std::cout << "a = " << a << ", b = " << b << std::endl;
//      g_c[0][c] += (2.*a*lambda_c[0][c]/(h_+1.e-6) - cp_*rho[0][c]*(2.*a*xc[2] + b)*(dhdt*xc[2] - B_w))/(h_+1.e-6) * cv[0][c];
//          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      g_c[0][c] += 2.*a*lambda_c[0][c]/(h_+1.e-6)/(h_+1.e-6) * cv[0][c];
      //          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      // manufactured solution without non-conservative part
//      g_c[0][c] -= (-2.*a*lambda_c[0][c]/(h_+1.e-6) -
//          cp_*rho[0][c]/(h_+1.e-6)*( dhdt*(3.*a*xc[2]*xc[2] + 2.*b*xc[2] + c) - B_w*(2.*a*xc[2] + b) ) ) * cv[0][c];
//          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      g_c[0][c] += ( 2.*a*lambda_c[0][c]/(h_+1.e-6)/(h_+1.e-6) - (2.*a*xc[2] + b) ) * cv[0][c];
//       g_c[0][c] += ( 2.*a*lambda_c[0][c]/(h_+1.e-6)/(h_+1.e-6) - cp_*rho[0][c]*dhdt/(h_+1.e-6)*(3.*a*xc[2]*xc[2] + 2.*b*xc[2] + c) + cp_*rho[0][c]*B_w/(h_+1.e-6)*(2.*a*xc[2] + b) ) * cv[0][c]
//                 + cp_*rho[0][c]*dhdt/(h_+1.e-6)*(a*xc[2]*xc[2] + b*xc[2] + c) * cv[0][c];
//      std::cout << "h_ = " << h_ << std::endl;
//      std::cout << "SRC = " << g_c[0][c] << std::endl;
       *
       */
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "Adding external source term" << std::endl;
    //    db_->WriteVector("  Q_ext", S->GetFieldData(source_key_).ptr(), false);
    db_->WriteVector("res (src)", g, false);
  }
}


void Lake_Thermo_PK::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
//  // external sources of energy (temperature dependent source)
//  if (is_source_term_ && is_source_term_differentiable_ &&
//      S->GetFieldEvaluator(source_key_)->IsDependency(S, key_)) {
//
//    Teuchos::RCP<const CompositeVector> dsource_dT;
//    if (!is_source_term_finite_differentiable_) {
//      // evaluate the derivative through the dag
//      S->GetFieldEvaluator(source_key_)->HasFieldDerivativeChanged(S, name_, key_);
//      dsource_dT = S->GetFieldData(Keys::getDerivKey(source_key_, key_));
//    } else {
//      // evaluate the derivative through finite differences
//      double eps = 1.e-8;
//      S->GetFieldData(key_, name_)->Shift(eps);
//      ChangedSolution();
//      S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
//      auto dsource_dT_nc = Teuchos::rcp(new CompositeVector(*S->GetFieldData(source_key_)));
//
//      S->GetFieldData(key_, name_)->Shift(-eps);
//      ChangedSolution();
//      S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
//
//      dsource_dT_nc->Update(-1/eps, *S->GetFieldData(source_key_), 1/eps);
//      dsource_dT = dsource_dT_nc;
//    }
//    db_->WriteVector("  dQ_ext/dT", dsource_dT.ptr(), false);
//    preconditioner_acc_->AddAccumulationTerm(*dsource_dT, -1.0, "cell", true);
//  }
}

// -------------------------------------------------------------
// Plug enthalpy into the boundary faces manually.
// This will be removed once boundary faces exist.
// -------------------------------------------------------------
void Lake_Thermo_PK::ApplyDirichletBCsToTemperature_(const Teuchos::Ptr<State>& S) {
  // put the boundary fluxes in faces for Dirichlet BCs.
  //  S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S, name_);

  const Epetra_MultiVector& temp_bf =
      *S->GetFieldData(temperature_key_)->ViewComponent("boundary_face",false);
  const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
  const Epetra_Map& face_map = mesh_->face_map(false);

  int nbfaces = temp_bf.MyLength();
  for (int bf=0; bf!=nbfaces; ++bf) {
    AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));

    if (bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
      bc_adv_->bc_value()[f] = temp_bf[0][bf];
    }
  }
}

void Lake_Thermo_PK::ApplyDirichletBCsToEnergy_(const Teuchos::Ptr<State>& S) {
  // put the boundary fluxes in faces for Dirichlet BCs.
  //  S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S, name_);

  S->GetFieldEvaluator(energy_key_)->HasFieldChanged(S, name_);
  const Epetra_MultiVector& enrg_bf =
      *S->GetFieldData(energy_key_)->ViewComponent("boundary_face",false);
  const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
  const Epetra_Map& face_map = mesh_->face_map(false);

  int nbfaces = enrg_bf.MyLength();
  for (int bf=0; bf!=nbfaces; ++bf) {
    AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));

    if (bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
      bc_adv_->bc_value()[f] = enrg_bf[0][bf];
    }
  }
}



} //namespace LakeThermo
} //namespace Amanzi
