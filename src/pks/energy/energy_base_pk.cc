/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

#include "energy_bc_factory.hh"
#include "advection_factory.hh"

#include "PDE_DiffusionFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_AdvectionUpwind.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"

#include "enthalpy_evaluator.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"
#include "pk_helpers.hh"
#include "Mesh_Algorithms.hh"

#include "energy_base.hh"

#define MORE_DEBUG_FLAG 0


namespace Amanzi {
namespace Energy {


EnergyBase::EnergyBase(Teuchos::ParameterList& FElist,
                       const Teuchos::RCP<Teuchos::ParameterList>& plist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution)
  : PK(FElist, plist, S, solution),
    PK_PhysicalBDF_Default(FElist, plist, S, solution),
    modify_predictor_with_consistent_faces_(false),
    modify_predictor_for_freezing_(false),
    coupled_to_subsurface_via_temp_(false),
    coupled_to_subsurface_via_flux_(false),
    coupled_to_surface_via_temp_(false),
    coupled_to_surface_via_flux_(false),
    decoupled_from_subsurface_(false),
    niter_(0),
    flux_exists_(true)
{
  // set a default error tolerance
  if (domain_.find("surface") != std::string::npos) {
    mass_atol_ = plist_->get<double>("mass absolute error tolerance", .01 * 55000.);
    soil_atol_ = 0.;
  } else {
    mass_atol_ = plist_->get<double>("mass absolute error tolerance", .5 * .1 * 55000.);
    soil_atol_ = 0.5 * 2000. * 620.e-6;
    // porosity * particle density soil * heat capacity soil * 1 degree
    // or, dry bulk density soil * heat capacity soil * 1 degree, in MJ
  }
  if (!plist_->isParameter("absolute error tolerance")) {
    plist_->set("absolute error tolerance", 76.e-6);
    // energy of 1 degree C of water per mass_atol, in MJ/mol water
  }

  // source terms
  is_source_term_ = plist_->get<bool>("source term", false);
  if (is_source_term_ && source_key_.empty()) {
    source_key_ = Keys::readKey(*plist_, domain_, "source", "total_energy_source");
  }
  is_source_term_differentiable_ = plist_->get<bool>("source term is differentiable", true);
  is_source_term_finite_differentiable_ = plist_->get<bool>("source term finite difference", false);

  // get keys
  conserved_key_ = Keys::readKey(*plist_, domain_, "conserved quantity", "energy");
  wc_key_ = Keys::readKey(*plist_, domain_, "water content", "water_content");
  enthalpy_key_ = Keys::readKey(*plist_, domain_, "enthalpy", "enthalpy");
  flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");
  energy_flux_key_ =
    Keys::readKey(*plist_, domain_, "diffusive energy flux", "diffusive_energy_flux");
  adv_energy_flux_key_ =
    Keys::readKey(*plist_, domain_, "advected energy flux", "advected_energy_flux");
  conductivity_key_ =
    Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");
  uw_conductivity_key_ =
    Keys::readKey(*plist_, domain_, "upwinded thermal conductivity", "upwind_thermal_conductivity");
  uf_key_ = Keys::readKey(*plist_, domain_, "unfrozen fraction", "unfrozen_fraction");

  if (S_->IsDeformableMesh(domain_))
    deform_key_ = Keys::readKey(*plist_, domain_, "deformation indicator", "base_porosity");
}


// -------------------------------------------------------------
// Setup
// -------------------------------------------------------------
void
EnergyBase::Setup()
{
  PK_PhysicalBDF_Default::Setup();
  SetupEnergy_();
  SetupPhysicalEvaluators_();
};


void
EnergyBase::SetupPhysicalEvaluators_()
{
  // Get data and evaluators needed by the PK
  // -- energy, energy evaluator, and energy derivative
  requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    conserved_key_, tag_next_, key_, tag_next_);

  // -- conductivity
  requireAtNext(conductivity_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}

// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Energy-like PKs.
// -------------------------------------------------------------
void
EnergyBase::SetupEnergy_()
{
  // Get data for special-case entities.
  requireAtNext(cell_vol_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<double>("atmospheric_pressure", Tags::DEFAULT);

  // require data for source terms
  if (is_source_term_) {
    requireAtNext(source_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    if (is_source_term_differentiable_ && (!is_source_term_finite_differentiable_)) {
      // NOTE, the following line is commented out because of the bug, amanzi/ats#167
      //&& S_->GetEvaluator(source_key_, tag_next_).IsDifferentiableWRT(*S_, key_, tag_next_)) {
      // require derivative of source
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
          source_key_, tag_next_, key_, tag_next_)
        .SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      // NOTE, remove SetMesh/AddComponent lines after fixing amanzi/ats#167.
      // The mesh should get set by the evaluator, but when
      // the evaluator isn't actually differentiable, it
      // doesn't get done.
    }
  }

  // is dynamic mesh?  If so, get a key for indicating when the mesh has changed.
  if (!deform_key_.empty()) S_->RequireEvaluator(deform_key_, tag_next_);

  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(mesh_, bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_diff_flux_ = bc_factory.CreateDiffusiveFlux();
  bc_flux_ = bc_factory.CreateTotalFlux();

  bc_adv_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));

  // -- nonlinear coefficient
  std::string method_name =
    plist_->get<std::string>("upwind conductivity method", "arithmetic mean");
  if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_, tag_next_));
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_, tag_next_));
  } else {
    Errors::Message message;
    message << "Energy PK has no upwinding method named: " << method_name;
    Exceptions::amanzi_throw(message);
  }

  std::string coef_location = upwinding_->CoefficientLocation();
  if (coef_location == "upwind: face") {
    S_->Require<CompositeVector, CompositeVectorSpace>(uw_conductivity_key_, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S_->Require<CompositeVector, CompositeVectorSpace>(uw_conductivity_key_, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  } else {
    Errors::Message message("Unknown upwind coefficient location in energy.");
    Exceptions::amanzi_throw(message);
  }
  S_->GetRecordW(uw_conductivity_key_, tag_next_, name_).set_io_vis(false);

  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.Create(mfd_plist, mesh_, bc_);
  matrix_diff_->SetTensorCoefficient(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // -- create the operators for the preconditioner
  //    diffusion
  // NOTE: Can this be a clone of the primary operator? --etc
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  if (!mfd_pc_plist.isParameter("discretization primary"))
    mfd_pc_plist.set("discretization primary",
                     mfd_plist.get<std::string>("discretization primary"));
  if (!mfd_pc_plist.isParameter("discretization secondary") &&
      mfd_plist.isParameter("discretization secondary"))
    mfd_pc_plist.set("discretization secondary",
                     mfd_plist.get<std::string>("discretization secondary"));
  if (!mfd_pc_plist.isParameter("schema") && mfd_plist.isParameter("schema"))
    mfd_pc_plist.set("schema", mfd_plist.get<Teuchos::Array<std::string>>("schema"));
  if (mfd_pc_plist.get<bool>("include Newton correction", false)) {
    if (mfd_pc_plist.get<std::string>("discretization primary") == "fv: default") {
      mfd_pc_plist.set("Newton correction", "true Jacobian");
    } else {
      mfd_pc_plist.set("Newton correction", "approximate Jacobian");
    }
  }

  precon_used_ = plist_->isSublist("preconditioner") || plist_->isSublist("inverse") ||
                 plist_->isSublist("linear solver");
  if (precon_used_) {
    auto& inv_list = mfd_pc_plist.sublist("inverse");
    inv_list.setParameters(plist_->sublist("inverse"));
    // old style... deprecate me!
    inv_list.setParameters(plist_->sublist("preconditioner"));
    inv_list.setParameters(plist_->sublist("linear solver"));
  }

  preconditioner_diff_ = opfactory.Create(mfd_pc_plist, mesh_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  //    If using approximate Jacobian for the preconditioner, we also
  //    need derivative information.  This means upwinding the
  //    derivative.
  jacobian_ = mfd_pc_plist.get<std::string>("Newton correction", "none") != "none";
  if (jacobian_) {
    // require the derivative dcond/dp
    S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        conductivity_key_, tag_next_, key_, tag_next_)
      .SetGhosted();

    if (mfd_pc_plist.get<std::string>("discretization primary") != "fv: default") {
      // MFD or NLFV -- upwind required
      duw_conductivity_key_ = Keys::getDerivKey(uw_conductivity_key_, key_);

      S_->Require<CompositeVector, CompositeVectorSpace>(duw_conductivity_key_, tag_next_, name_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

      upwinding_deriv_ =
        Teuchos::rcp(new Operators::UpwindTotalFlux(name_, tag_next_, energy_flux_key_, 1.e-8));
    } else {
      // FV -- no upwinding
      duw_conductivity_key_ = "";
    }
  } else {
    duw_conductivity_key_ = "";
  }

  // -- accumulation terms
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ =
    Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  //  -- advection terms
  // -- set up the evaluator for enthalpy, whether or not we advect the thing,
  //    we may need it for exchange fluxes with surface/subsurface coupling,
  //    etc.
  if (plist_->isSublist("enthalpy evaluator")) {
    Teuchos::ParameterList& enth_list = S_->GetEvaluatorList(enthalpy_key_);
    enth_list.setParameters(plist_->sublist("enthalpy evaluator"));
    enth_list.set<std::string>("evaluator type", "enthalpy");
  }
  is_advection_term_ = plist_->get<bool>("include thermal advection", true);
  if (is_advection_term_) {
    // -- create the forward operator for the advection term
    Teuchos::ParameterList advect_plist = plist_->sublist("advection");
    matrix_adv_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(advect_plist, mesh_));
    matrix_adv_->SetBCs(bc_adv_, bc_adv_);

    implicit_advection_ = !plist_->get<bool>("explicit advection", false);
    if (implicit_advection_) {
      // -- operator for preconditioner
      implicit_advection_in_pc_ =
        !plist_->get<bool>("supress advective terms in preconditioner", false);
      if (implicit_advection_in_pc_) {
        Teuchos::ParameterList advect_plist_pc = plist_->sublist("advection preconditioner");
        preconditioner_adv_ =
          Teuchos::rcp(new Operators::PDE_AdvectionUpwind(advect_plist_pc, preconditioner_));
        preconditioner_adv_->SetBCs(bc_adv_, bc_adv_);
      }

      // -- enthalpy data, evaluator
      requireAtNext(enthalpy_key_, tag_next_, *S_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
        ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
          enthalpy_key_, tag_next_, key_, tag_next_)
        .SetGhosted();

      // -- diagnostic for advected energy
      requireAtNext(adv_energy_flux_key_, tag_next_, *S_, name_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

    } else {
      // -- enthalpy data, evaluator
      requireAtCurrent(enthalpy_key_, tag_current_, *S_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
        ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

      // -- diagnostic for advected energy
      requireAtCurrent(adv_energy_flux_key_, tag_current_, *S_, name_)
        .SetMesh(mesh_)
        ->SetGhosted()
        ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    }
  }

  // coupling terms
  // -- coupled to a surface via a Neumann condition
  coupled_to_surface_via_flux_ = plist_->get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    std::string domain_surf;
    if (domain_ == "domain" || domain_ == "") {
      domain_surf = plist_->get<std::string>("surface domain name", "surface");
    } else {
      domain_surf = plist_->get<std::string>("surface domain name", "surface_" + domain_);
    }
    ss_flux_key_ = Keys::readKey(
      *plist_, domain_surf, "surface-subsurface energy flux", "surface_subsurface_energy_flux");
    S_->Require<CompositeVector, CompositeVectorSpace>(ss_flux_key_, tag_next_)
      .SetMesh(S_->GetMesh(domain_surf))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- coupled to a surface via a Dirichlet condition
  coupled_to_surface_via_temp_ = plist_->get<bool>("coupled to surface via temperature", false);
  if (coupled_to_surface_via_temp_) {
    std::string domain_surf;
    if (domain_ == "domain" || domain_ == "") {
      domain_surf = plist_->get<std::string>("surface domain name", "surface");
    } else {
      domain_surf = plist_->get<std::string>("surface domain name", "surface_" + domain_);
    }
    ss_primary_key_ = Keys::readKey(*plist_, domain_surf, "temperature", "temperature");
    S_->Require<CompositeVector, CompositeVectorSpace>(ss_primary_key_, tag_next_)
      .SetMesh(S_->GetMesh(domain_surf))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  decoupled_from_subsurface_ =
    plist_->get<bool>("decoupled from subsurface", false); //surface-only system

  // -- Make sure coupling isn't flagged multiple ways.
  if (coupled_to_surface_via_flux_ && coupled_to_surface_via_temp_) {
    Errors::Message message(
      "Energy PK requested both flux and temperature coupling -- choose one.");
    Exceptions::amanzi_throw(message);
  }

  // Require the primary variable
  CompositeVectorSpace matrix_cvs = matrix_->RangeMap();
  matrix_cvs.AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_)
    .Update(matrix_cvs)
    ->SetGhosted();

  // require a water flux field
  requireAtNext(flux_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(flux_key_, tag_current_, *S_, name_);

  // require a water content field -- used for computing energy density in the
  // error norm
  requireAtNext(wc_key_, tag_next_, *S_).SetMesh(mesh_)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtCurrent(wc_key_, tag_current_, *S_, name_);

  // Require fields for the energy fluxes
  requireAtNext(energy_flux_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // Globalization and other timestep control flags
  modify_predictor_for_freezing_ = plist_->get<bool>("modify predictor for freezing", false);
  modify_predictor_with_consistent_faces_ =
    plist_->get<bool>("modify predictor with consistent faces", false);
  T_limit_ = plist_->get<double>("limit correction to temperature change [K]", -1.);
};


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void
EnergyBase::Initialize()
{
  // initialize BDF stuff and physical domain stuff
  PK_PhysicalBDF_Default::Initialize();

  // initialize energy fluxes
  S_->GetW<CompositeVector>(energy_flux_key_, tag_next_, name()).PutScalar(0.0);
  S_->GetRecordW(energy_flux_key_, tag_next_, name()).set_initialized();
  changedEvaluatorPrimary(energy_flux_key_, tag_next_, *S_);

  if (is_advection_term_) {
    Tag adv_energy_flux_tag = implicit_advection_ ? tag_next_ : tag_current_;
    S_->GetW<CompositeVector>(adv_energy_flux_key_, adv_energy_flux_tag, name()).PutScalar(0.0);
    S_->GetRecordW(adv_energy_flux_key_, adv_energy_flux_tag, name()).set_initialized();
    changedEvaluatorPrimary(adv_energy_flux_key_, adv_energy_flux_tag, *S_);
  }

  // initialize upwind conductivities
  S_->GetW<CompositeVector>(uw_conductivity_key_, tag_next_, name()).PutScalar(0.0);
  S_->GetRecordW(uw_conductivity_key_, tag_next_, name()).set_initialized();
  if (!duw_conductivity_key_.empty()) {
    S_->GetW<CompositeVector>(duw_conductivity_key_, tag_next_, name()).PutScalar(0.);
    S_->GetRecordW(duw_conductivity_key_, tag_next_, name()).set_initialized();
  }
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void
EnergyBase::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  // saves primary variable
  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, tag_next);

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;

  assign(wc_key_, tag_current, tag_next, *S_);
  assign(flux_key_, tag_current, tag_next, *S_);
}


bool
EnergyBase::UpdateConductivityData_(const Tag& tag)
{
  Teuchos::RCP<const CompositeVector> cond = S_->GetPtr<CompositeVector>(conductivity_key_, tag);
  bool update = S_->GetEvaluator(conductivity_key_, tag).Update(*S_, name_);

  if (update) {
    Teuchos::RCP<CompositeVector> uw_cond =
      S_->GetPtrW<CompositeVector>(uw_conductivity_key_, tag, name_);
    upwinding_->Update(*cond, *uw_cond, *S_);

    if (uw_cond->HasComponent("face")) uw_cond->ScatterMasterToGhosted("face");
  }
  return update;
}


bool
EnergyBase::UpdateConductivityDerivativeData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating conductivity derivatives? ";

  bool update = S_->GetEvaluator(conductivity_key_, tag).UpdateDerivative(*S_, name_, key_, tag);

  if (update) {
    Teuchos::RCP<const CompositeVector> dcond =
      S_->GetDerivativePtr<CompositeVector>(conductivity_key_, tag, key_, tag);
    if (!duw_conductivity_key_.empty()) {
      // get upwind conductivity data
      Teuchos::RCP<CompositeVector> duw_cond =
        S_->GetPtrW<CompositeVector>(duw_conductivity_key_, tag, name_);
      duw_cond->PutScalar(0.);
      upwinding_deriv_->Update(*dcond, *duw_cond, *S_);
      if (duw_cond->HasComponent("face")) duw_cond->ScatterMasterToGhosted("face");
    } else {
      dcond->ScatterMasterToGhosted("cell");
    }
  }
  return update;
}

// -----------------------------------------------------------------------------
// Compute boundary condition functionals at the tag
// -----------------------------------------------------------------------------
void
EnergyBase::ComputeBoundaryConditions_(const Tag& tag)
{
  bc_temperature_->Compute(S_->get_time(tag));
  bc_diff_flux_->Compute(S_->get_time(tag));
  bc_flux_->Compute(S_->get_time(tag));
}

// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void
EnergyBase::UpdateBoundaryConditions_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating BCs." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  auto& adv_markers = bc_adv_->bc_model();
  auto& adv_values = bc_adv_->bc_value();

  for (unsigned int n = 0; n != markers.size(); ++n) {
    markers[n] = Operators::OPERATOR_BC_NONE;
    values[n] = 0.0;
    adv_markers[n] = Operators::OPERATOR_BC_NONE;
    adv_values[n] = 0.0;
  }

  // Dirichlet temperature boundary conditions
  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc != bc_temperature_->end();
       ++bc) {
    AmanziMesh::Entity_ID_List cells;
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    values[f] = bc->second;
    adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  }

  // Neumann flux boundary conditions
  for (Functions::BoundaryFunction::Iterator bc = bc_flux_->begin(); bc != bc_flux_->end(); ++bc) {
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_NEUMANN;
    values[f] = bc->second;
    adv_markers[f] = Operators::OPERATOR_BC_NEUMANN;
    // push all onto diffusion, assuming that the incoming enthalpy is 0 (likely water flux is 0)
  }

  // Neumann diffusive flux, not Neumann TOTAL flux.  Potentially advective flux.
  for (Functions::BoundaryFunction::Iterator bc = bc_diff_flux_->begin();
       bc != bc_diff_flux_->end();
       ++bc) {
    int f = bc->first;
    markers[f] = Operators::OPERATOR_BC_NEUMANN;
    values[f] = bc->second;
    adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  }

  // Dirichlet temperature boundary conditions from a coupled surface.
  if (coupled_to_surface_via_temp_) {
    // Face is Dirichlet with value of surface temp
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh(Keys::getDomain(ss_flux_key_));
    const Epetra_MultiVector& temp =
      *S_->GetPtr<CompositeVector>("surface_temperature", tag)->ViewComponent("cell", false);

    int ncells_surface = temp.MyLength();
    for (int c = 0; c != ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f = surface->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

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
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh(Keys::getDomain(ss_flux_key_));
    const Epetra_MultiVector& flux =
      *S_->Get<CompositeVector>(ss_flux_key_, tag).ViewComponent("cell", false);

    int ncells_surface = flux.MyLength();
    for (int c = 0; c != ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f = surface->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

      // -- set that value to Neumann
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      // flux provided by the coupler is in units of J / s, whereas Neumann BCs are J/s/A
      values[f] = flux[0][c] / mesh_->getFaceArea(f);

      // -- mark advective BCs as Dirichlet: this ensures the surface
      //    temperature is picked up and advective fluxes are treated
      //    via advection operator, not diffusion operator.
      adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }

  // mark all remaining boundary conditions as zero diffusive flux conditions
  int nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (markers[f] == Operators::OPERATOR_BC_NONE) {
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      if (ncells == 1) {
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
        adv_markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      }
    }
  }

  // set the face temperature on boundary faces
  auto temp = S_->GetPtrW<CompositeVector>(key_, tag, name_);
  applyDirichletBCs(*bc_, *temp);
};


// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool
EnergyBase::IsAdmissible(Teuchos::RCP<const TreeVector> up)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Checking admissibility..." << std::endl;

  // For some reason, wandering PKs break most frequently with an unreasonable
  // temperature.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> temp = up->Data();
  double minT, maxT;

  const Epetra_MultiVector& temp_c = *temp->ViewComponent("cell", false);
  double minT_c(1.e6), maxT_c(-1.e6);
  int min_c(-1), max_c(-1);
  for (int c = 0; c != temp_c.MyLength(); ++c) {
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
    const Epetra_MultiVector& temp_f = *temp->ViewComponent("face", false);
    for (int f = 0; f != temp_f.MyLength(); ++f) {
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
  mesh_->getComm()->MaxAll(&maxT_l, &maxT, 1);
  mesh_->getComm()->MinAll(&minT_l, &minT, 1);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "    Admissible T? (min/max): " << minT << ",  " << maxT << std::endl;
  }

  Teuchos::RCP<const Comm_type> comm_p = mesh_->getComm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();

  if (minT < 200.0 || maxT > 330.0) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << " is not admissible, as it is not within bounds of constitutive models:"
                 << std::endl;
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
        const Epetra_MultiVector& temp_f = *temp->ViewComponent("face", false);
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
bool
EnergyBase::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Modifying predictor:" << std::endl;

  // update boundary conditions
  //ComputeBoundaryConditions_(tag_next_);
  bc_temperature_->Compute(S_->get_time(tag_next_));
  bc_flux_->Compute(S_->get_time(tag_next_));
  UpdateBoundaryConditions_(tag_next_);

  // predictor modification
  bool modified = false;
  if (modify_predictor_for_freezing_) {
    const Epetra_MultiVector& u0_c = *u0->Data()->ViewComponent("cell", false);
    Epetra_MultiVector& u_c = *u->Data()->ViewComponent("cell", false);

    for (int c = 0; c != u0_c.MyLength(); ++c) {
      if (u0_c[0][c] > 273.15 && u_c[0][c] < 273.15) {
        u_c[0][c] = 273.15 - .00001;
        modified = true;
      } else if (u0_c[0][c] < 273.15 && u_c[0][c] > 273.15) {
        u_c[0][c] = 273.15 + .00001;
        modified = true;
      }
    }

    if (u0->Data()->HasComponent("boundary_face")) {
      const Epetra_MultiVector& u0_bf = *u0->Data()->ViewComponent("boundary_face", false);
      Epetra_MultiVector& u_bf = *u->Data()->ViewComponent("boundary_face", false);

      for (int bf = 0; bf != u0_bf.MyLength(); ++bf) {
        if (u0_bf[0][bf] > 273.15 && u_bf[0][bf] < 273.15) {
          u_bf[0][bf] = 273.15 - .00001;
          modified = true;
        } else if (u0_bf[0][bf] < 273.15 && u_bf[0][bf] > 273.15) {
          u_bf[0][bf] = 273.15 + .00001;
          modified = true;
        }
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
void
EnergyBase::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u)
{
  if (!u->HasComponent("face")) return;

  // average cells to faces to give a reasonable initial guess
  u->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell", true);
  Epetra_MultiVector& u_f = *u->ViewComponent("face", false);

  int f_owned = u_f.MyLength();
  for (int f = 0; f != f_owned; ++f) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n = 0; n != ncells; ++n) { face_value += u_c[0][cells[n]]; }
    u_f[0][f] = face_value / ncells;
  }
  ChangedSolution();

  // use old BCs
  // // update boundary conditions
  // bc_temperature_->Compute(S_->get_time(tag_next_));
  // bc_diff_flux_->Compute(S_->get_time(tag_next_));
  // bc_flux_->Compute(S_->get_time(tag_next_));
  // UpdateBoundaryConditions_(tag_next_);

  // use old conductivity
  // div K_e grad u
  //  bool update = UpdateConductivityData_(tag_next_);
  Teuchos::RCP<const CompositeVector> conductivity =
    S_->GetPtr<CompositeVector>(uw_conductivity_key_, tag_next_);

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
EnergyBase::ModifyCorrection(double h,
                             Teuchos::RCP<const TreeVector> res,
                             Teuchos::RCP<const TreeVector> u,
                             Teuchos::RCP<TreeVector> du)
{
  // update diffusive flux correction -- note this is not really modifying the correction as far as NKA is concerned
  const auto& T_vec = *u->Data();
  if (T_vec.HasComponent("boundary_face")) {
    const Epetra_MultiVector& T_bf = *T_vec.ViewComponent("boundary_face", false);
    const Epetra_MultiVector& T_c = *T_vec.ViewComponent("cell", false);

    const Epetra_MultiVector& dT_c = *du->Data()->ViewComponent("cell", false);
    Epetra_MultiVector& dT_bf = *du->Data()->ViewComponent("boundary_face", false);

    for (int bf = 0; bf != T_bf.MyLength(); ++bf) {
      AmanziMesh::Entity_ID f = getBoundaryFaceFace(*mesh_, bf);

      // NOTE: this should get refactored into a helper class, much like predictor_delegate_bc_flux
      // as this would be necessary to deal with general discretizations.  Note that this is not
      // needed in cases where boundary faces are already up to date (e.g. MFD, maybe NLFV?)
      if (bc_markers()[f] == Operators::OPERATOR_BC_NEUMANN &&
          bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
        // diffusive flux BC
        AmanziMesh::Entity_ID c = getFaceOnBoundaryInternalCell(*mesh_, f);
        const auto& Acc = matrix_diff_->local_op()->matrices_shadow[f];
        double T_bf_val =
          (Acc(0, 0) * (T_c[0][c] - dT_c[0][c]) - bc_values()[f] * mesh_->getFaceArea(f)) / Acc(0, 0);
        dT_bf[0][bf] = T_bf[0][bf] - T_bf_val;
      }
    }
  }


  // max temperature correction
  int my_limited = 0;
  int n_limited = 0;
  if (T_limit_ > 0.) {
    for (CompositeVector::name_iterator comp = du->Data()->begin(); comp != du->Data()->end();
         ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp, false);

      double max;
      du_c.NormInf(&max);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Max temperature correction (" << *comp << ") = " << max << std::endl;
      }

      for (int c = 0; c != du_c.MyLength(); ++c) {
        if (std::abs(du_c[0][c]) > T_limit_) {
          du_c[0][c] = ((du_c[0][c] > 0) - (du_c[0][c] < 0)) * T_limit_;
          my_limited++;
        }
      }
    }
    mesh_->getComm()->SumAll(&my_limited, &n_limited, 1);
  }

  if (n_limited > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) { *vo_->os() << "  limited by temperature." << std::endl; }
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  }
  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

} // namespace Energy
} // namespace Amanzi
