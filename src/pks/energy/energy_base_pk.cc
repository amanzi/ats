/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */
#include "boost/algorithm/string/predicate.hpp"

#include "energy_bc_factory.hh"
#include "advection_factory.hh"

#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_AdvectionUpwind.hh"
#include "LinearOperatorFactory.hh"
#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"
#include "enthalpy_evaluator.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "energy_base.hh"

#define MORE_DEBUG_FLAG 0


namespace Amanzi {
namespace Energy {


EnergyBase::EnergyBase(Teuchos::ParameterList& FElist,
                       const Teuchos::RCP<Teuchos::ParameterList>& plist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    PK_PhysicalBDF_Default(FElist, plist, S, solution),
    modify_predictor_with_consistent_faces_(false),
    modify_predictor_for_freezing_(false),
    coupled_to_subsurface_via_temp_(false),
    coupled_to_subsurface_via_flux_(false),
    coupled_to_surface_via_temp_(false),
    coupled_to_surface_via_flux_(false),
    niter_(0),
    flux_exists_(true),
    implicit_advection_(true) {

  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "energy");

  // set a default error tolerance
  if (domain_.find("surface") != std::string::npos) {
    mass_atol_ = plist_->get<double>("mass absolute error tolerance",
                                     .01 * 55000.);
    soil_atol_ = 0.;
  } else {
    mass_atol_ = plist_->get<double>("mass absolute error tolerance",
                                     .5 * .1 * 55000.);
    soil_atol_ = 0.5 * 2000. * 620.e-6;  // porosity * particle density soil * heat capacity soil * 1 degree
                    // or, dry bulk density soil * heat capacity soil * 1 degree, in MJ
  }

  if (!plist_->isParameter("absolute error tolerance")) {
    plist_->set("absolute error tolerance", 76.e-6); // energy of 1 degree C of water per mass_atol, in MJ/mol water
  }
}


// -------------------------------------------------------------
// Setup
// -------------------------------------------------------------
void EnergyBase::Setup(const Teuchos::Ptr<State>& S) {
  PK_PhysicalBDF_Default::Setup(S);

  SetupEnergy_(S);
  SetupPhysicalEvaluators_(S);

};


void EnergyBase::SetupEnergy_(const Teuchos::Ptr<State>& S) {
  // Set up keys if they were not already set.
  energy_key_ = Keys::readKey(*plist_, domain_, "energy", "energy");
  wc_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");
  enthalpy_key_ = Keys::readKey(*plist_, domain_, "enthalpy", "enthalpy");
  flux_key_ = Keys::readKey(*plist_, domain_, "mass flux", "mass_flux");
  energy_flux_key_ = Keys::readKey(*plist_, domain_, "diffusive energy flux", "diffusive_energy_flux");
  adv_energy_flux_key_ = Keys::readKey(*plist_, domain_, "advected energy flux", "advected_energy_flux");
  conductivity_key_ = Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");
  uw_conductivity_key_ = Keys::readKey(*plist_, domain_, "upwinded thermal conductivity", "upwind_thermal_conductivity");

  // Get data for special-case entities.
  S->RequireField(cell_vol_key_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(cell_vol_key_);
  S->RequireScalar("atmospheric_pressure");

  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(mesh_, bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_diff_flux_ = bc_factory.CreateDiffusiveFlux();
  bc_flux_ = bc_factory.CreateTotalFlux();

  bc_adv_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // -- nonlinear coefficient
  std::string method_name = plist_->get<std::string>("upwind conductivity method",
          "arithmetic mean");
  if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            conductivity_key_, uw_conductivity_key_));
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            conductivity_key_, uw_conductivity_key_));
  } else {
    std::stringstream messagestream;
    messagestream << "Energy PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  std::string coef_location = upwinding_->CoefficientLocation();
  if (coef_location == "upwind: face") {  
    S->RequireField(uw_conductivity_key_, name_)->SetMesh(mesh_)
        ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S->RequireField(uw_conductivity_key_, name_)->SetMesh(mesh_)
        ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  } else {
    Errors::Message message("Unknown upwind coefficient location in energy.");
    Exceptions::amanzi_throw(message);
  }
  S->GetField(uw_conductivity_key_,name_)->set_io_vis(false);
  
  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.Create(mfd_plist, mesh_, bc_);
  matrix_diff_->SetTensorCoefficient(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // -- create the forward operator for the advection term
  Teuchos::ParameterList advect_plist = plist_->sublist("advection");
  matrix_adv_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(advect_plist, mesh_));
  matrix_adv_->SetBCs(bc_adv_, bc_adv_);
  
  // -- create the operators for the preconditioner
  //    diffusion
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  if (!mfd_pc_plist.isParameter("discretization primary"))
    mfd_pc_plist.set("discretization primary", mfd_plist.get<std::string>("discretization primary"));
  if (!mfd_pc_plist.isParameter("discretization secondary") && mfd_plist.isParameter("discretization secondary"))
    mfd_pc_plist.set("discretization secondary", mfd_plist.get<std::string>("discretization secondary"));
  if (!mfd_pc_plist.isParameter("schema") && mfd_plist.isParameter("schema"))
    mfd_pc_plist.set("schema", mfd_plist.get<Teuchos::Array<std::string> >("schema"));

  preconditioner_diff_ = opfactory.Create(mfd_pc_plist, mesh_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();
  
  //    If using approximate Jacobian for the preconditioner, we also
  //    need derivative information.  This means upwinding the
  //    derivative.
  jacobian_ = mfd_pc_plist.get<std::string>("Newton correction", "none") != "none";
  if (jacobian_) {
    // if (preconditioner_->RangeMap().HasComponent("face")) {
    if (mfd_pc_plist.get<std::string>("discretization primary") != "fv: default"){
      // MFD -- upwind required
      dconductivity_key_ = Keys::getDerivKey(conductivity_key_, key_);
      duw_conductivity_key_ = Keys::getDerivKey(uw_conductivity_key_, key_);
        
      S->RequireField(duw_conductivity_key_, name_)
        ->SetMesh(mesh_)->SetGhosted()
        ->SetComponent("face", AmanziMesh::FACE, 1);

      // upwinding_deriv_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
      //                                 dconductivity_key_, duw_conductivity_key_));
      upwinding_deriv_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                                      dconductivity_key_, duw_conductivity_key_,
                                      energy_flux_key_, 1.e-8));

    } else {
      // FV -- no upwinding
      dconductivity_key_ = Keys::getDerivKey(conductivity_key_, key_);
      duw_conductivity_key_ = "";
    }
  } else {
    dconductivity_key_ = "";
    duw_conductivity_key_ = "";
  }
  

  // -- accumulation terms
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  //  -- advection terms
  implicit_advection_ = !plist_->get<bool>("explicit advection", false);
  if (implicit_advection_) {
    implicit_advection_in_pc_ = !plist_->get<bool>("supress advective terms in preconditioner", false);

    if (implicit_advection_in_pc_) {
      Teuchos::ParameterList advect_plist_pc = plist_->sublist("advection preconditioner");
      preconditioner_adv_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(advect_plist_pc, preconditioner_));
      preconditioner_adv_->SetBCs(bc_adv_, bc_adv_);
    }
  }

  //    symbolic assemble
  precon_used_ = plist_->isSublist("preconditioner");
  if (precon_used_) {
    preconditioner_->SymbolicAssembleMatrix();
    preconditioner_->InitializePreconditioner(plist_->sublist("preconditioner"));

    //    Potentially create a linear solver
    if (plist_->isSublist("linear solver")) {
      Teuchos::ParameterList linsolve_sublist = plist_->sublist("linear solver");
      AmanziSolvers::LinearOperatorFactory<Operators::Operator,CompositeVector,CompositeVectorSpace> fac;
      lin_solver_ = fac.Create(linsolve_sublist, preconditioner_);
    } else {
      lin_solver_ = preconditioner_;
    }
  }  

  // -- advection of enthalpy
  S->RequireField(enthalpy_key_)->SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  Teuchos::ParameterList enth_plist = plist_->sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", enthalpy_key_);
  Teuchos::RCP<EnthalpyEvaluator> enth =
    Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S->SetFieldEvaluator(enthalpy_key_, enth);

  // source terms
  is_source_term_ = plist_->get<bool>("source term");
  is_source_term_differentiable_ = plist_->get<bool>("source term is differentiable", true);
  is_source_term_finite_differentiable_ = plist_->get<bool>("source term finite difference", false);
  if (is_source_term_) {
    source_key_ = Keys::readKey(*plist_, domain_, "source", "total_energy_source");
    S->RequireField(source_key_)->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(source_key_);
  }

  // coupling terms
  // -- subsurface PK, coupled to the surface
  coupled_to_surface_via_flux_ =
      plist_->get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    std::string domain_surf;
    if (domain_ == "domain" || domain_ == "") {
      domain_surf = plist_->get<std::string>("surface domain name", "surface");
    } else {
      domain_surf = plist_->get<std::string>("surface domain name", "surface_"+domain_);
    }
    ss_flux_key_ = Keys::readKey(*plist_, domain_surf, "surface-subsurface energy flux", "surface_subsurface_energy_flux");
    S->RequireField(ss_flux_key_)
        ->SetMesh(S->GetMesh(domain_surf))
        ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  coupled_to_surface_via_temp_ =
      plist_->get<bool>("coupled to surface via temperature", false);
  if (coupled_to_surface_via_temp_) {
    // surface temperature used for BCs
    S->RequireField("surface_temperature")
        ->SetMesh(S->GetMesh("surface"))
        ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- Make sure coupling isn't flagged multiple ways.
  if (coupled_to_surface_via_flux_ && coupled_to_surface_via_temp_) {
    Errors::Message message("Energy PK requested both flux and temperature coupling -- choose one.");
    Exceptions::amanzi_throw(message);
  }

  CompositeVectorSpace matrix_cvs = matrix_->RangeMap();
  matrix_cvs.AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1); 

  // -- primary variable
  S->RequireField(key_, name_)->Update(matrix_cvs)->SetGhosted();
  
  // Require a field for the mass flux for advection.
  flux_exists_ = S->HasField(flux_key_); // this bool is needed to know if PK
                                         // makes flux or we need an
                                         // independent variable evaluator

  S->RequireField(flux_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("face", AmanziMesh::FACE, 1);

  // Require a field for the energy fluxes.
  S->RequireField(energy_flux_key_, name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField(adv_energy_flux_key_, name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);

  // -- simply limit to close to 0
  modify_predictor_for_freezing_ =
      plist_->get<bool>("modify predictor for freezing", false);
  // -- ewc and other predictors can result in odd face values
  modify_predictor_with_consistent_faces_ =
      plist_->get<bool>("modify predictor with consistent faces", false);

  // correction controls
  T_limit_ = plist_->get<double>("limit correction to temperature change [K]", -1.);
};


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void EnergyBase::Initialize(const Teuchos::Ptr<State>& S) {
  // initialize BDF stuff and physical domain stuff
  PK_PhysicalBDF_Default::Initialize(S);

#if MORE_DEBUG_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << i;
    S->GetFieldData(namestream.str(),name_)->PutScalar(0.);
    S->GetField(namestream.str(),name_)->set_initialized();

    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << i;
    S->GetFieldData(solnstream.str(),name_)->PutScalar(0.);
    S->GetField(solnstream.str(),name_)->set_initialized();
  }

#endif

  // initialize energy flux
  S->GetFieldData(energy_flux_key_, name_)->PutScalar(0.0);
  S->GetField(energy_flux_key_, name_)->set_initialized();
  S->GetFieldData(adv_energy_flux_key_, name_)->PutScalar(0.0);
  S->GetField(adv_energy_flux_key_, name_)->set_initialized();
  S->GetFieldData(uw_conductivity_key_, name_)->PutScalar(0.0);
  S->GetField(uw_conductivity_key_, name_)->set_initialized();
  if (!duw_conductivity_key_.empty()) {
    S->GetFieldData(duw_conductivity_key_, name_)->PutScalar(0.);
    S->GetField(duw_conductivity_key_, name_)->set_initialized();
  }
  
  // potentially initialize mass flux
  if (!flux_exists_) {
    S->GetField(flux_key_, name_)->Initialize(plist_->sublist(flux_key_));
  }

};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void EnergyBase::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {

  double dt = t_new - t_old;
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;
  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, S);

  bc_temperature_->Compute(S->time());
  bc_diff_flux_->Compute(S->time());
  bc_flux_->Compute(S->time());
  UpdateBoundaryConditions_(S.ptr());
  
  niter_ = 0;
  bool update = UpdateConductivityData_(S.ptr());
  update |= S->GetFieldEvaluator(key_)->HasFieldChanged(S.ptr(), name_);

  if (update) {
    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);
    Teuchos::RCP<const CompositeVector> conductivity = S->GetFieldData(uw_conductivity_key_);
    matrix_diff_->global_operator()->Init();
    matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
    matrix_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
    matrix_diff_->ApplyBCs(true, true, true);

    Teuchos::RCP<CompositeVector> eflux = S->GetFieldData(energy_flux_key_, name_);
    matrix_diff_->UpdateFlux(temp.ptr(), eflux.ptr());

    // calculate the advected energy as a diagnostic
    Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_key_);
    matrix_adv_->Setup(*flux);
    S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S.ptr(), name_);
    Teuchos::RCP<const CompositeVector> enth = S->GetFieldData(enthalpy_key_);;
    ApplyDirichletBCsToEnthalpy_(S.ptr());

    Teuchos::RCP<CompositeVector> adv_energy = S->GetFieldData(adv_energy_flux_key_, name_);  
    matrix_adv_->UpdateFlux(enth.ptr(), flux.ptr(), bc_adv_, adv_energy.ptr());
  }
};

/* ******************************************************************
* Returns the first cell attached to a boundary face.   
****************************************************************** */
int EnergyBase::BoundaryFaceGetCell(int f) const
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  return cells[0];
}

/* ******************************************************************
* For FV methods, set the boundary face temp.
****************************************************************** */
void EnergyBase::ApplyDirichletBCsToBoundaryFace_(const Teuchos::Ptr<CompositeVector>& temp)
{

  Epetra_MultiVector& temp_bf = *temp->ViewComponent("boundary_face", false);
  const Epetra_MultiVector& temp_c = *temp->ViewComponent("cell", false);
  const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
  const Epetra_Map& face_map = mesh_->face_map(false);

  const std::vector<int>& bc_model = bc_->bc_model();
  const std::vector<double>& bc_value = bc_->bc_value();
  
  int nbfaces = temp_bf.MyLength();

  for (int bf=0; bf!=nbfaces; ++bf) {
    AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      temp_bf[0][bf] = bc_value[f];
    } else {
      temp_bf[0][bf] = temp_c[0][BoundaryFaceGetCell(f)];
    }
  }
}


bool EnergyBase::UpdateConductivityData_(const Teuchos::Ptr<State>& S) {
  bool update = S->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S, name_);
  if (update) {
    upwinding_->Update(S);

    Teuchos::RCP<CompositeVector> uw_cond =
        S->GetFieldData(uw_conductivity_key_, name_);
    if (uw_cond->HasComponent("face"))
      uw_cond->ScatterMasterToGhosted("face");
  }
  return update;
}


bool EnergyBase::UpdateConductivityDerivativeData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating conductivity derivatives? ";

  bool update = S->GetFieldEvaluator(conductivity_key_)
    ->HasFieldDerivativeChanged(S, name_, key_);

  if (update) {
    if (!duw_conductivity_key_.empty()) {
      upwinding_deriv_->Update(S);

      Teuchos::RCP<CompositeVector> duw_cond =
        S->GetFieldData(duw_conductivity_key_, name_);
      if (duw_cond->HasComponent("face"))
        duw_cond->ScatterMasterToGhosted("face");
    } else {
      Teuchos::RCP<const CompositeVector> dcond =
        S->GetFieldData(dconductivity_key_);
      dcond->ScatterMasterToGhosted("cell");
    }
  }
  return update;
}

// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void EnergyBase::UpdateBoundaryConditions_(
    const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  auto& adv_markers = bc_adv_->bc_model();
  auto& adv_values = bc_adv_->bc_value();
  
  for (unsigned int n=0; n!=markers.size(); ++n) {
    markers[n] = Operators::OPERATOR_BC_NONE;
    values[n] = 0.0;
    adv_markers[n] = Operators::OPERATOR_BC_NONE;
    adv_values[n] = 0.0;
  }

  // Dirichlet temperature boundary conditions
  for (Functions::BoundaryFunction::Iterator bc=bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    AmanziMesh::Entity_ID_List cells;
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    values[f] = bc->second;
    adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  }

  // Neumann flux boundary conditions
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_NEUMANN;
    values[f] = bc->second;
    adv_markers[f] = Operators::OPERATOR_BC_NEUMANN;
    // push all onto diffusion, assuming that the incoming enthalpy is 0 (likely mass flux is 0)
  }

  // Neumann diffusive flux, not Neumann TOTAL flux.  Potentially advective flux.
  for (Functions::BoundaryFunction::Iterator bc=bc_diff_flux_->begin();
       bc!=bc_diff_flux_->end(); ++bc) {
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_NEUMANN;
    values[f] = bc->second;
    adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  }
  
  // Dirichlet temperature boundary conditions from a coupled surface.
  if (coupled_to_surface_via_temp_) {
    // Face is Dirichlet with value of surface temp
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh(Keys::getDomain(ss_flux_key_));
    const Epetra_MultiVector& temp = *S->GetFieldData("surface_temperature")
        ->ViewComponent("cell",false);

    int ncells_surface = temp.MyLength();
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to dirichlet
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = temp[0][c];
      adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }

  // surface coupling
  if (coupled_to_surface_via_flux_) {
    // Diffusive fluxes are given by the residual of the surface equation.
    // Advective fluxes are given by the surface temperature and whatever flux we have.
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh(Keys::getDomain(ss_flux_key_));
    const Epetra_MultiVector& flux =
      *S->GetFieldData(ss_flux_key_)->ViewComponent("cell",false);

    int ncells_surface = flux.MyLength();
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to Neumann
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      // flux provided by the coupler is in units of J / s, whereas Neumann BCs are J/s/A
      values[f] = flux[0][c] / mesh_->face_area(f);

      // -- mark advective BCs as Dirichlet: this ensures the surface
      //    temperature is picked up and advective fluxes are treated
      //    via advection operator, not diffusion operator.
      adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }

  // mark all remaining boundary conditions as zero diffusive flux conditions
  AmanziMesh::Entity_ID_List cells;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (markers[f] == Operators::OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
        adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      }
    }
  }

  // set the face temperature on boundary faces
  auto temp = S->GetFieldData(key_, name_);
  if (temp->HasComponent("boundary_face")) {
    ApplyDirichletBCsToBoundaryFace_(temp.ptr());
  }
};



// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool EnergyBase::IsAdmissible(Teuchos::RCP<const TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Checking admissibility..." << std::endl;

  // For some reason, wandering PKs break most frequently with an unreasonable
  // temperature.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> temp = up->Data();
  double minT, maxT;
  
  const Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
  double minT_c(1.e6), maxT_c(-1.e6);
  int min_c(-1), max_c(-1);
  for (int c=0; c!=temp_c.MyLength(); ++c) {
    if (temp_c[0][c] < minT_c) {
      minT_c = temp_c[0][c];
      min_c = c;
    }
    if (temp_c[0][c] > maxT_c) {
      maxT_c = temp_c[0][c];
      max_c = c;
    }
  }

  double minT_f(1.e6), maxT_f(-1.e6);
  int min_f(-1), max_f(-1);
  if (temp->HasComponent("face")) {
    const Epetra_MultiVector& temp_f = *temp->ViewComponent("face",false);
    for (int f=0; f!=temp_f.MyLength(); ++f) {
      if (temp_f[0][f] < minT_f) {
        minT_f = temp_f[0][f];
        min_f = f;
      }
      if (temp_f[0][f] > maxT_f) {
        maxT_f = temp_f[0][f];
        max_f = f;
      }
    }
    minT = std::min(minT_c, minT_f);
    maxT = std::max(maxT_c, maxT_f);

  } else {
    minT = minT_c;
    maxT = maxT_c;
  }

  double minT_l = minT;
  double maxT_l = maxT;
  mesh_->get_comm()->MaxAll(&maxT_l, &maxT, 1);
  mesh_->get_comm()->MinAll(&minT_l, &minT, 1);
  
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "    Admissible T? (min/max): " << minT << ",  " << maxT << std::endl;
  }

  Teuchos::RCP<const Comm_type> comm_p = mesh_->get_comm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();


  
  if (minT < 200.0 || maxT > 300.0) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << " is not admissible, as it is not within bounds of constitutive models:" << std::endl;
      ENorm_t global_minT_c, local_minT_c;
      ENorm_t global_maxT_c, local_maxT_c;

      local_minT_c.value = minT_c;
      local_minT_c.gid = temp_c.Map().GID(min_c);
      local_maxT_c.value = maxT_c;
      local_maxT_c.gid = temp_c.Map().GID(max_c);

      MPI_Allreduce(&local_minT_c, &global_minT_c, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
      MPI_Allreduce(&local_maxT_c, &global_maxT_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

      *vo_->os() << "   cells (min/max): [" << global_minT_c.gid << "] " << global_minT_c.value
                 << ", [" << global_maxT_c.gid << "] " << global_maxT_c.value << std::endl;

      if (temp->HasComponent("face")) {
        const Epetra_MultiVector& temp_f = *temp->ViewComponent("face",false);
        ENorm_t global_minT_f, local_minT_f;
        ENorm_t global_maxT_f, local_maxT_f;

        local_minT_f.value = minT_f;
        local_minT_f.gid = temp_f.Map().GID(min_f);
        local_maxT_f.value = maxT_f;
        local_maxT_f.gid = temp_f.Map().GID(max_f);
        

        MPI_Allreduce(&local_minT_f, &global_minT_f, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
        MPI_Allreduce(&local_maxT_f, &global_maxT_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
        *vo_->os() << "   cells (min/max): [" << global_minT_f.gid << "] " << global_minT_f.value
                   << ", [" << global_maxT_f.gid << "] " << global_maxT_f.value << std::endl;
      }
    }
    return false;
  }
  return true;
}


// -----------------------------------------------------------------------------
// BDF takes a prediction step -- make sure it is physical and otherwise ok.
// -----------------------------------------------------------------------------
bool EnergyBase::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_.ptr());
  
  // push Dirichlet data into predictor
  if (u->Data()->HasComponent("boundary_cell")) {
    ApplyBoundaryConditions_(u->Data().ptr());
  }

  bool modified = false;
  if (modify_predictor_for_freezing_) {
    const Epetra_MultiVector& u0_c = *u0->Data()->ViewComponent("cell",false);
    Epetra_MultiVector& u_c = *u->Data()->ViewComponent("cell",false);

    for (int c=0; c!=u0_c.MyLength(); ++c) {
      if (u0_c[0][c] > 273.15 && u_c[0][c] < 273.15) {
        u_c[0][c] = 273.15 - .00001;
        modified = true;
      }
    }
  }

  if (modify_predictor_with_consistent_faces_) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  modifications for consistent face temperatures." << std::endl;
    CalculateConsistentFaces(u->Data().ptr());
    modified = true;
  }
  return modified;
}


// -----------------------------------------------------------------------------
// Given an arbitrary set of cell values, calculate consitent face constraints.
//
//  This is useful for prediction steps, hacky preconditioners, etc.
// -----------------------------------------------------------------------------
void EnergyBase::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {

  // average cells to faces to give a reasonable initial guess
  u->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell",true);
  Epetra_MultiVector& u_f = *u->ViewComponent("face",false);

  int f_owned = u_f.MyLength();
  for (int f=0; f!=f_owned; ++f) {
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += u_c[0][cells[n]];
    }
    u_f[0][f] = face_value / ncells;
  }
  ChangedSolution();
  
  // use old BCs
  // // update boundary conditions
  // bc_temperature_->Compute(S_next_->time());
  // bc_diff_flux_->Compute(S_next_->time());
  // bc_flux_->Compute(S_next_->time());
  // UpdateBoundaryConditions_(S_next_.ptr());

  // use old conductivity
  // div K_e grad u
  //  bool update = UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  // Update the preconditioner
  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, u);
  //matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  matrix_diff_->ApplyBCs(true, true, true);

  // derive the consistent faces, involves a solve
  matrix_diff_->UpdateConsistentFaces(*u);
}



AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
EnergyBase::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                             Teuchos::RCP<const TreeVector> u,
                             Teuchos::RCP<TreeVector> du) {

  int my_limited = 0;
  int n_limited = 0;
  if (T_limit_ > 0.) {
    for (CompositeVector::name_iterator comp=du->Data()->begin();
         comp!=du->Data()->end(); ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);

      double max;
      du_c.NormInf(&max);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Max temperature correction (" << *comp << ") = " << max << std::endl;
      }

      for (int c=0; c!=du_c.MyLength(); ++c) {
        if (std::abs(du_c[0][c]) > T_limit_) {
          du_c[0][c] = ((du_c[0][c] > 0) - (du_c[0][c] < 0)) * T_limit_;
          my_limited++;
        }
      }
    }
    mesh_->get_comm()->SumAll(&my_limited, &n_limited, 1);
  }

  if (n_limited > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "  limited by temperature." << std::endl;
    }
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  }
  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

} // namespace Energy
} // namespace Amanzi
