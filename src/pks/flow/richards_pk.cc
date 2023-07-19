/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Epetra_Import.h"

#include "flow_bc_factory.hh"

#include "Point.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "predictor_delegate_bc_flux.hh"
#include "wrm_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content_evaluator.hh"
#include "OperatorDefs.hh"
#include "BoundaryFlux.hh"
#include "pk_helpers.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------

Richards::Richards(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, glist, S, solution),
    PK_PhysicalBDF_Default(pk_tree, glist, S, solution),
    coupled_to_surface_via_head_(false),
    coupled_to_surface_via_flux_(false),
    infiltrate_only_if_unfrozen_(false),
    modify_predictor_with_consistent_faces_(false),
    modify_predictor_wc_(false),
    modify_predictor_bc_flux_(false),
    modify_predictor_first_bc_flux_(false),
    upwind_from_prev_flux_(false),
    clobber_boundary_flux_dir_(false),
    vapor_diffusion_(false),
    perm_scale_(1.),
    jacobian_(false),
    jacobian_lag_(0),
    iter_(0),
    iter_counter_time_(0.),
    fixed_kr_(false)
{
  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .5 * .1 * 55000.); // phi * s * nl

  // get field names
  conserved_key_ = Keys::readKey(*plist_, domain_, "conserved", "water_content");
  mass_dens_key_ = Keys::readKey(*plist_, domain_, "mass density", "mass_density_liquid");
  molar_dens_key_ = Keys::readKey(*plist_, domain_, "molar density", "molar_density_liquid");
  perm_key_ = Keys::readKey(*plist_, domain_, "permeability", "permeability");
  coef_key_ = Keys::readKey(*plist_, domain_, "conductivity", "relative_permeability");
  uw_coef_key_ =
    Keys::readKey(*plist_, domain_, "upwinded conductivity", "upwind_relative_permeability");
  flux_key_ = Keys::readKey(*plist_, domain_, "darcy flux", "water_flux");
  flux_dir_key_ = Keys::readKey(*plist_, domain_, "darcy flux direction", "water_flux_direction");
  velocity_key_ = Keys::readKey(*plist_, domain_, "darcy velocity", "darcy_velocity");
  sat_key_ = Keys::readKey(*plist_, domain_, "saturation", "saturation_liquid");
  sat_gas_key_ = Keys::readKey(*plist_, domain_, "saturation gas", "saturation_gas");
  sat_ice_key_ = Keys::readKey(*plist_, domain_, "saturation ice", "saturation_ice");
  capillary_pressure_gas_liq_key_ =
    Keys::readKey(*plist_, domain_, "capillary_pressure_gas_liq", "capillary_pressure_gas_liq");
  capillary_pressure_liq_ice_key_ =
    Keys::readKey(*plist_, domain_, "capillary_pressure_liq_ice", "capillary_pressure_liq_ice");

  if (S_->IsDeformableMesh(domain_))
    deform_key_ = Keys::readKey(*plist_, domain_, "deformation indicator", "base_porosity");

  // all manipulation of evaluator lists should happen in constructors (pre-setup)
  // -- WRM: This deals with deprecated location for the WRM list (in the PK).
  if (plist_->isSublist("water retention evaluator")) {
    auto& wrm_plist = S_->GetEvaluatorList(sat_key_);
    wrm_plist.setParameters(plist_->sublist("water retention evaluator"));
    wrm_plist.set("evaluator type", "WRM");
  }
  if (S_->GetEvaluatorList(coef_key_).numParams() == 0) {
    Teuchos::ParameterList& kr_plist = S_->GetEvaluatorList(coef_key_);
    kr_plist.setParameters(S_->GetEvaluatorList(sat_key_));
    kr_plist.set<std::string>("evaluator type", "WRM rel perm");
  }

  // scaling for permeability for better "nondimensionalization"
  perm_scale_ = plist_->get<double>("permeability rescaling", 1.e7);
  S_->GetEvaluatorList(coef_key_).set<double>("permeability rescaling", perm_scale_);
}

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void
Richards::Setup()
{
  PK_PhysicalBDF_Default::Setup();
  SetupRichardsFlow_();
  SetupPhysicalEvaluators_();
};


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Richards-like PKs.
// -------------------------------------------------------------
void
Richards::SetupRichardsFlow_()
{
  // Get data for special-case entities.
  S_->Require<CompositeVector, CompositeVectorSpace>(cell_vol_key_, tag_next_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(cell_vol_key_, tag_next_);

  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_head_ = bc_factory.CreateHead();
  bc_level_ = bc_factory.CreateFixedLevel();
  bc_flux_ = bc_factory.CreateMassFlux();
  bc_seepage_ = bc_factory.CreateSeepageFacePressure();
  bc_seepage_->Compute(0.); // compute at t=0 to set up
  std::tie(bc_seepage_infilt_explicit_, bc_seepage_infilt_) =
    bc_factory.CreateSeepageFacePressureWithInfiltration();
  bc_seepage_infilt_->Compute(0.); // compute at t=0 to set up
  bc_rho_water_ = bc_plist.get<double>("hydrostatic water density [kg m^-3]", 1000.);

  // -- linear tensor coefficients
  // permeability type - scalar or tensor?
  Teuchos::ParameterList& perm_list = S_->GetEvaluatorList(perm_key_);
  std::string perm_type = perm_list.get<std::string>("permeability type", "scalar");
  if (perm_type == "scalar") {
    perm_tensor_rank_ = 1;
    num_perm_vals_ = 1;
  } else if (perm_type == "horizontal and vertical") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = 2;
  } else if (perm_type == "diagonal tensor") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = mesh_->getSpaceDimension();
  } else if (perm_type == "full tensor") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = (mesh_->getSpaceDimension() == 3) ? 6 : 3;
  } else {
    Errors::Message message(
      "`permeability type` must be one of the following: \"scalar\", \"diagonal tensor\", \"full "
      "tensor\", or \"horizontal and vertical\".");
    Exceptions::amanzi_throw(message);
  }

  // is dynamic mesh?  If so, get a key for indicating when the mesh has changed.
  if (!deform_key_.empty()) S_->RequireEvaluator(deform_key_, tag_next_);

  // data allocation -- move to State!
  unsigned int c_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (unsigned int c = 0; c != c_owned; ++c) {
    (*K_)[c].Init(mesh_->getSpaceDimension(), perm_tensor_rank_);
  }

  // -- nonlinear coefficients/upwinding
  // if coupled to the surface, how do we deal with the surface face
  if (plist_->isParameter("surface rel perm strategy")) {
    clobber_policy_ = plist_->get<std::string>("surface rel perm strategy");
  } else if (plist_->get<bool>("clobber surface rel perm", false)) {
    clobber_policy_ = "clobber";
  } else if (plist_->get<bool>("max surface rel perm", false)) {
    clobber_policy_ = "max";
  } else if (plist_->get<bool>("unsaturated clobber surface rel perm", false)) {
    clobber_policy_ = "unsaturated";
  } else {
    clobber_policy_ = "none";
  }
  clobber_boundary_flux_dir_ =
    plist_->get<bool>("clobber boundary flux direction for upwinding", false);

  // what upwinding method to use
  std::string method_name =
    plist_->get<std::string>("relative permeability method", "upwind with Darcy flux");
  if (method_name == "upwind with gravity") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_, tag_next_, K_));
    Krel_method_ = Operators::UPWIND_METHOD_GRAVITY;
  } else if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_, tag_next_));
    Krel_method_ = Operators::UPWIND_METHOD_CENTERED;
  } else if (method_name == "upwind with Darcy flux") {
    double flux_eps = plist_->get<double>("upwind flux epsilon", 1.e-5);
    upwinding_ =
      Teuchos::rcp(new Operators::UpwindTotalFlux(name_, tag_next_, flux_dir_key_, flux_eps));
    Krel_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_, tag_next_));
    Krel_method_ = Operators::UPWIND_METHOD_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Richards Flow PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // require the data on appropriate locations
  std::string coef_location = upwinding_->CoefficientLocation();
  RequireNonlinearCoefficient_(uw_coef_key_, coef_location);

  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  mfd_plist.set("gravity", true);

  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.CreateWithGravity(mfd_plist, mesh_, bc_);
  matrix_ = matrix_diff_->global_operator();

  // -- create the operator, data for flux directions
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = opfactory.CreateWithGravity(face_diff_list, mesh_, bc_);

  S_->Require<CompositeVector, CompositeVectorSpace>(flux_dir_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // -- create the operators for the preconditioner
  //    diffusion
  // NOTE: Can this be a clone of the primary operator? --etc
  //      get the discretiation type
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  mfd_pc_plist.set("gravity", true);
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

  //    get the inverse method
  precon_used_ = plist_->isSublist("preconditioner") || plist_->isSublist("inverse") ||
                 plist_->isSublist("linear solver");
  if (precon_used_) {
    mfd_pc_plist.set("inverse", plist_->sublist("inverse"));
    // old style... deprecate me!
    mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("preconditioner"));
    mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("linear solver"));
  }

  //    create the operator
  preconditioner_diff_ = opfactory.CreateWithGravity(mfd_pc_plist, mesh_, bc_);
  preconditioner_ = preconditioner_diff_->global_operator();

  //    If using approximate Jacobian for the preconditioner, we also need
  //    derivative information.  For now this means upwinding the derivative.
  jacobian_ = mfd_pc_plist.get<std::string>("Newton correction", "none") != "none";
  if (jacobian_) {
    jacobian_lag_ = mfd_pc_plist.get<int>("Newton correction lag", 0);

    // require the derivative drel_perm/dp
    S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        coef_key_, tag_next_, key_, tag_next_)
      .SetGhosted();
    if (mfd_pc_plist.get<std::string>("discretization primary") != "fv: default") {
      // MFD -- upwind required, require data
      duw_coef_key_ = Keys::getDerivKey(uw_coef_key_, key_);

      // note this is always done on faces using total flux upwinding
      RequireNonlinearCoefficient_(duw_coef_key_, "upwind: face");

      // note, this is here to be consistent -- unclear whether the 1.e-3 is useful or not?
      double flux_eps = plist_->get<double>("upwind flux epsilon", 1.e-5);
      upwinding_deriv_ = Teuchos::rcp(
        new Operators::UpwindTotalFlux(name_, tag_next_, flux_dir_key_, 1.e-3 * flux_eps));
    } else {
      // FV -- no upwinding of derivative
      duw_coef_key_ = std::string();
    }
  }

  // -- accumulation terms
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set<std::string>("entity kind", "cell");
  preconditioner_acc_ =
    Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  // // -- vapor diffusion terms
  // vapor_diffusion_ = plist_->get<bool>("include vapor diffusion", false);
  // if (vapor_diffusion_){
  //   AMANZI_ASSERT(0); // untested!

  //   // Create the vapor diffusion vectors
  //   S_->Require<CompositeVector,CompositeVectorSpace>("vapor_diffusion_pressure", tag_next_,  name_).SetMesh(mesh_)->SetGhosted()
  //       ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  //   S_->GetRecordW("vapor_diffusion_pressure",name_)->set_io_vis(true);

  //   S_->Require<CompositeVector,CompositeVectorSpace>("vapor_diffusion_temperature", tag_next_,  name_).SetMesh(mesh_)->SetGhosted()
  //     ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  //   S_->GetRecordW("vapor_diffusion_temperature",name_)->set_io_vis(true);

  //   // operator for the vapor diffusion terms
  //   matrix_vapor_ = Operators::CreateMatrixMFD(mfd_plist, mesh_);
  // }

  // -- source terms
  is_source_term_ = plist_->get<bool>("source term", false);
  if (is_source_term_) {
    if (source_key_.empty())
      source_key_ = Keys::readKey(*plist_, domain_, "source", "water_source");
    source_term_is_differentiable_ = plist_->get<bool>("source term is differentiable", true);
    explicit_source_ = plist_->get<bool>("explicit source term", false);

    requireAtNext(source_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    if (source_term_is_differentiable_) {
      // require derivative of source
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        source_key_, tag_next_, key_, tag_next_);
    }
  }

  // coupling to the surface
  // -- coupling done by a Neumann condition
  coupled_to_surface_via_flux_ = plist_->get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    Key domain_surf = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");
    ss_flux_key_ =
      Keys::readKey(*plist_, domain_surf, "surface-subsurface flux", "surface_subsurface_flux");
    S_->Require<CompositeVector, CompositeVectorSpace>(ss_flux_key_, tag_next_)
      .SetMesh(S_->GetMesh(domain_surf))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- coupling done by a Dirichlet condition
  coupled_to_surface_via_head_ = plist_->get<bool>("coupled to surface via head", false);
  if (coupled_to_surface_via_head_) {
    Key domain_surf = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");
    ss_primary_key_ = Keys::readKey(*plist_, domain_surf, "pressure", "pressure");
    S_->Require<CompositeVector, CompositeVectorSpace>(ss_primary_key_, tag_next_)
      .SetMesh(S_->GetMesh(domain_surf))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- Make sure coupling isn't flagged multiple ways.
  if (coupled_to_surface_via_flux_ && coupled_to_surface_via_head_) {
    Errors::Message message("Richards PK requested both flux and head coupling -- choose one.");
    Exceptions::amanzi_throw(message);
  }

  // Require fields and evaluators
  // -- pressure, the primary variable
  //  NOTE: no need to require evaluator for p here, either at the old or new
  //  times, as this was done in pk_physical.  All we have to do is set the
  //  structure.
  CompositeVectorSpace matrix_cvs = matrix_->RangeMap();
  compute_boundary_values_ = plist_->get<bool>("compute boundary values", false);
  if (compute_boundary_values_)
    matrix_cvs.AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_)
    .Update(matrix_cvs)
    ->SetGhosted();

  // -- flux is managed here as a primary variable
  requireAtNext(flux_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // -- also need a velocity, but only for vis/diagnostics, so might as well
  // -- only keep at NEXT
  requireAtNext(velocity_key_, Tags::NEXT, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 3);

  // Globalization and other timestep control flags
  // -- predictors
  modify_predictor_with_consistent_faces_ =
    plist_->get<bool>("modify predictor with consistent faces", false);
  modify_predictor_bc_flux_ = plist_->get<bool>("modify predictor for flux BCs", false);
  modify_predictor_first_bc_flux_ =
    plist_->get<bool>("modify predictor for initial flux BCs", false);
  modify_predictor_wc_ = plist_->get<bool>("modify predictor via water content", false);

  // -- correctors
  p_limit_ = plist_->get<double>("limit correction to pressure change [Pa]", -1.);
  patm_limit_ =
    plist_->get<double>("limit correction to pressure change when crossing atmospheric [Pa]", -1.);

  // -- valid step controls
  sat_change_limit_ = plist_->get<double>("max valid change in saturation in a time step [-]", -1.);
  sat_ice_change_limit_ =
    plist_->get<double>("max valid change in ice saturation in a time step [-]", -1.);
}

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void
Richards::SetupPhysicalEvaluators_()
{
  // -- Absolute permeability.
  requireAtNext(perm_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_perm_vals_);

  // -- water content, and evaluator, and derivative for PC
  requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    conserved_key_, tag_next_, key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(conserved_key_, tag_current_, *S_, name_);

  // -- Water retention evaluators
  // -- saturation
  requireAtNext(sat_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireAtNext(sat_gas_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  auto& wrm = S_->RequireEvaluator(sat_key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(sat_key_, tag_current_, *S_, name_);

  // -- rel perm
  requireAtNext(coef_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

  // -- get the WRM models
  auto wrm_eval = dynamic_cast<Flow::WRMEvaluator*>(&wrm);
  AMANZI_ASSERT(wrm_eval != nullptr);
  wrms_ = wrm_eval->get_WRMs();

  // -- molar density used to infer liquid Darcy velocity from flux
  requireAtNext(molar_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // -- liquid mass density for the gravity fluxes
  requireAtNext(mass_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


// -------------------------------------------------------------
// Helper function and customization point for upwinded coefs.
// -------------------------------------------------------------
void
Richards::RequireNonlinearCoefficient_(const Key& key, const std::string& coef_location)
{
  if (coef_location == "upwind: face") {
    S_->Require<CompositeVector, CompositeVectorSpace>(key, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S_->Require<CompositeVector, CompositeVectorSpace>(key, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  } else {
    Errors::Message message("Unknown upwind coefficient location in Richards flow.");
    Exceptions::amanzi_throw(message);
  }
  S_->GetRecordW(key, tag_next_, name_).set_io_vis(false);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void
Richards::Initialize()
{
  // Initialize via hydrostatic balance
  if (!S_->GetRecordW(key_, tag_next_, name_).initialized()) InitializeHydrostatic_(tag_next_);

  // Initialize in the standard ways
  PK_PhysicalBDF_Default::Initialize();

  // Set extra fields as initialized -- these don't currently have evaluators,
  // and will be initialized in the call to commit_state()
  S_->GetW<CompositeVector>(uw_coef_key_, tag_next_, name_).PutScalar(1.0);
  S_->GetRecordW(uw_coef_key_, tag_next_, name_).set_initialized();

  if (!duw_coef_key_.empty()) {
    S_->GetW<CompositeVector>(duw_coef_key_, tag_next_, name_).PutScalar(1.0);
    S_->GetRecordW(duw_coef_key_, tag_next_, name_).set_initialized();
  }

  S_->GetW<CompositeVector>(flux_key_, tag_next_, name()).PutScalar(0.0);
  S_->GetRecordW(flux_key_, tag_next_, name()).set_initialized();
  changedEvaluatorPrimary(flux_key_, tag_next_, *S_);

  S_->GetW<CompositeVector>(flux_dir_key_, tag_next_, name()).PutScalar(0.0);
  S_->GetRecordW(flux_dir_key_, tag_next_, name()).set_initialized();
  S_->GetW<CompositeVector>(velocity_key_, Tags::NEXT, name()).PutScalar(0.0);
  S_->GetRecordW(velocity_key_, Tags::NEXT, name()).set_initialized();

  // absolute perm
  SetAbsolutePermeabilityTensor_(tag_next_);

  // operators
  const AmanziGeometry::Point& g = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  matrix_diff_->SetGravity(g);
  matrix_diff_->SetBCs(bc_, bc_);
  matrix_diff_->SetTensorCoefficient(K_);

  preconditioner_diff_->SetGravity(g);
  preconditioner_diff_->SetBCs(bc_, bc_);
  preconditioner_diff_->SetTensorCoefficient(K_);

  face_matrix_diff_->SetGravity(g);
  face_matrix_diff_->SetBCs(bc_, bc_);
  face_matrix_diff_->SetTensorCoefficient(K_);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);

  // if (vapor_diffusion_){
  //   //vapor diffusion
  //   matrix_vapor_->CreateMFDmassMatrices(Teuchos::null);
  //   // residual vector for vapor diffusion
  //   res_vapor = Teuchos::rcp(new CompositeVector(*S_->GetPtr<CompositeVector>(key_)));
  // }
};


void
Richards::InitializeHydrostatic_(const Tag& tag)
{
  // constant head over the surface
  if (plist_->sublist("initial condition").isParameter("hydrostatic head [m]")) {
    double head_wt = plist_->sublist("initial condition").get<double>("hydrostatic head [m]");
    double rho =
      plist_->sublist("initial condition").get<double>("hydrostatic water density [kg m^-3]");
    int ncols = mesh_->columns.num_columns_owned;

    const auto& gvec = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    int z_index = mesh_->getSpaceDimension() - 1;
    double g = -gvec[z_index];
    double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    // set pressure on the column of faces and cells
    Teuchos::RCP<CompositeVector> pres = S_->GetPtrW<CompositeVector>(key_, tag, name());
    Teuchos::RCP<Epetra_IntVector> flags = Teuchos::null;
    if (pres->HasComponent("face"))
      flags = Teuchos::rcp(new Epetra_IntVector(*pres->Map().Map("face", false)));

    { // context for viewcomponent -- do cells
      Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
      Teuchos::RCP<Epetra_MultiVector> pres_f = Teuchos::null;
      if (pres->HasComponent("face")) { pres_f = pres->ViewComponent("face", false); }

      int ncells_per = -1;
      for (int col = 0; col != ncols; ++col) {
        auto col_cells = mesh_->columns.getCells(col);
        if (ncells_per < 0) ncells_per = col_cells.size();
        AMANZI_ASSERT(col_cells.size() == ncells_per);
        auto col_faces = mesh_->columns.getFaces(col);
        AMANZI_ASSERT(col_faces.size() == col_cells.size() + 1);
        double z_wt = mesh_->getFaceCentroid(col_faces[0])[z_index] + head_wt;

        if (pres_f.get()) {
          (*pres_f)[0][col_faces[0]] = p_atm + rho * g * head_wt;
          (*flags)[col_faces[0]] = 1;
        }
        for (int lcv_c = 0; lcv_c != col_cells.size(); ++lcv_c) {
          AmanziMesh::Entity_ID c = col_cells[lcv_c];
          AmanziMesh::Entity_ID f = col_faces[lcv_c + 1];
          pres_c[0][c] = p_atm + rho * g * (z_wt - mesh_->getCellCentroid(c)[z_index]);
          if (pres_f.get()) {
            (*pres_f)[0][f] = p_atm + rho * g * (z_wt - mesh_->getFaceCentroid(f)[z_index]);
            (*flags)[f] = 1;
          }
        }
      }
    }

    if (pres->HasComponent("face")) {
      // communicate, then deal with horizontal-normal faces
      pres->ScatterMasterToGhosted("cell");
      {
        const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
        Epetra_MultiVector& pres_f = *pres->ViewComponent("face", false);
        for (AmanziMesh::Entity_ID f = 0; f != pres_f.MyLength(); ++f) {
          if (!(*flags)[f]) {
            auto f_cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
            if (f_cells.size() == 1) {
              // boundary face, use the cell value as the water table is
              // assumed to parallel the cell structure
              pres_f[0][f] = pres_c[0][f_cells[0]];
              (*flags)[f] = 1;
            } else {
              AMANZI_ASSERT(f_cells.size() == 2);
              // interpolate between cells
              pres_f[0][f] = (pres_c[0][f_cells[0]] + pres_c[0][f_cells[1]]) / 2.;
              (*flags)[f] = 1;
            }
          }
        }
      }
    }
    S_->GetRecordW(key_, tag, name()).set_initialized();
  }

  // constant head datum
  if (plist_->sublist("initial condition").isParameter("hydrostatic water level [m]")) {
    double z_wt = plist_->sublist("initial condition").get<double>("hydrostatic water level [m]");
    double rho =
      plist_->sublist("initial condition").get<double>("hydrostatic water density [kg m^-3]");

    int z_index = mesh_->getSpaceDimension() - 1;
    const auto& gravity = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    double g = -gravity[z_index];

    double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    Teuchos::RCP<CompositeVector> pres = S_->GetPtrW<CompositeVector>(key_, tag, name());
    Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
    for (int c = 0; c != pres_c.MyLength(); ++c) {
      pres_c[0][c] = p_atm + rho * g * (z_wt - mesh_->getCellCentroid(c)[z_index]);
    }

    if (pres->HasComponent("face")) {
      Epetra_MultiVector& pres_f = *pres->ViewComponent("face", false);
      for (int f = 0; f != pres_f.MyLength(); ++f) {
        pres_f[0][f] = p_atm + rho * g * (z_wt - mesh_->getFaceCentroid(f)[z_index]);
      }
    }
    S_->GetRecordW(key_, tag, name()).set_initialized();
  }
}


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void
Richards::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  // saves primary variable
  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, tag_next);

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;

  // also save saturation
  assign(sat_key_, tag_current, tag_next, *S_);
  if (S_->HasRecordSet(sat_ice_key_)) { assign(sat_ice_key_, tag_current, tag_next, *S_); }
};


// -----------------------------------------------------------------------------
// Check for controls on saturation
// -----------------------------------------------------------------------------
bool
Richards::ValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Validating time step." << std::endl;

  if (sat_change_limit_ > 0.0) {
    const Epetra_MultiVector& sl_new =
      *S_->GetPtr<CompositeVector>(sat_key_, tag_next_)->ViewComponent("cell", false);
    const Epetra_MultiVector& sl_old =
      *S_->GetPtr<CompositeVector>(sat_key_, tag_current_)->ViewComponent("cell", false);
    Epetra_MultiVector dsl(sl_new);
    dsl.Update(-1., sl_old, 1.);
    auto change = maxValLoc(*dsl(0));

    if (change.value > sat_change_limit_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max sl change=" << change.value
                   << " > limit=" << sat_change_limit_ << " at cell GID " << change.gid
                   << std::endl;
      return false;
    }
  }

  if (S_->HasRecordSet(sat_ice_key_) && (sat_ice_change_limit_ > 0.0)) {
    const Epetra_MultiVector& si_new =
      *S_->GetPtr<CompositeVector>(sat_ice_key_, tag_next_)->ViewComponent("cell", false);
    const Epetra_MultiVector& si_old =
      *S_->GetPtr<CompositeVector>(sat_ice_key_, tag_current_)->ViewComponent("cell", false);
    Epetra_MultiVector dsi(si_new);
    dsi.Update(-1., si_old, 1.);
    auto change = maxValLoc(*dsi(0));

    if (change.value > sat_ice_change_limit_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max si change=" << change.value
                   << " > limit=" << sat_ice_change_limit_ << " at cell GID " << change.gid
                   << std::endl;
      return false;
    }
  }
  return PK_PhysicalBDF_Default::ValidStep();
}


// -----------------------------------------------------------------------------
// Update any diagnostic variables prior to vis (in this case velocity field).
// -----------------------------------------------------------------------------
void
Richards::CalculateDiagnostics(const Tag& tag)
{
  AMANZI_ASSERT(tag == Tags::NEXT); // what else would this be?
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Calculating diagnostic variables." << std::endl;

  // update the cell velocities
  UpdateBoundaryConditions_(tag_next_);

  Teuchos::RCP<const CompositeVector> pres = S_->GetPtr<CompositeVector>(key_, tag_next_);
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_);
  Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag_next_);
  // update the stiffness matrix
  matrix_diff_->SetDensity(rho);
  matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // derive fluxes
  Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(flux_key_, tag_next_, name_);
  matrix_diff_->UpdateFlux(pres.ptr(), flux.ptr());
  UpdateVelocity_(tag);
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool
Richards::UpdatePermeabilityData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability?";
  if (fixed_kr_) return false;

  Teuchos::RCP<const CompositeVector> rel_perm = S_->GetPtr<CompositeVector>(coef_key_, tag);
  bool update_perm = S_->GetEvaluator(coef_key_, tag).Update(*S_, name_);

  // requirements due to the upwinding method
  if (Krel_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    bool update_dir = S_->GetEvaluator(mass_dens_key_, tag).Update(*S_, name_);
    update_dir |= S_->GetEvaluator(key_, tag).Update(*S_, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag);
      Teuchos::RCP<CompositeVector> flux_dir =
        S_->GetPtrW<CompositeVector>(flux_dir_key_, tag, name_);
      Teuchos::RCP<const CompositeVector> pres = S_->GetPtr<CompositeVector>(key_, tag);

      if (!deform_key_.empty() &&
          S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " flux dir"))
        face_matrix_diff_->SetTensorCoefficient(K_);
      face_matrix_diff_->SetDensity(rho);
      face_matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
      face_matrix_diff_->ApplyBCs(true, true, true);
      face_matrix_diff_->UpdateFlux(pres.ptr(), flux_dir.ptr());

      if (clobber_boundary_flux_dir_) {
        Epetra_MultiVector& flux_dir_f = *flux_dir->ViewComponent("face", false);

        auto& markers = bc_markers();
        auto& values = bc_values();

        for (int f = 0; f != markers.size(); ++f) {
          if (markers[f] == Operators::OPERATOR_BC_NEUMANN) {
            auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
            AMANZI_ASSERT(cells.size() == 1);
            int c = cells[0];
            const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
            int i = std::find(faces.begin(), faces.end(), f) - faces.begin();

            flux_dir_f[0][f] = values[f] * dirs[i];
          }
        }
      }
    }

    update_perm |= update_dir;
  }

  if (update_perm) {
    Teuchos::RCP<CompositeVector> uw_rel_perm =
      S_->GetPtrW<CompositeVector>(uw_coef_key_, tag, name_);

    // Move rel perm on boundary_faces into uw_rel_perm on faces
    const Epetra_Import& vandelay = mesh_->getBoundaryFaceImporter();
    const Epetra_MultiVector& rel_perm_bf = *rel_perm->ViewComponent("boundary_face", false);
    {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
      uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    }

    // Upwind, only overwriting boundary faces if the wind says to do so.
    upwinding_->Update(*rel_perm, *uw_rel_perm, *S_);

    if (clobber_policy_ == "clobber") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
      uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    } else if (clobber_policy_ == "max") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
      const auto& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
      const auto& bfmap = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
      for (int bf = 0; bf != rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        if (rel_perm_bf[0][bf] > uw_rel_perm_f[0][f]) { uw_rel_perm_f[0][f] = rel_perm_bf[0][bf]; }
      }
    } else if (clobber_policy_ == "unsaturated") {
      // clobber only when the interior cell is unsaturated
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
      const Epetra_MultiVector& pres =
        *S_->Get<CompositeVector>(key_, tag).ViewComponent("cell", false);
      const auto& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
      const auto& bfmap = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
      for (int bf = 0; bf != rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        auto fcells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        AMANZI_ASSERT(fcells.size() == 1);
        if (pres[0][fcells[0]] < 101225.) {
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf];
        } else if (pres[0][fcells[0]] < 101325.) {
          double frac = (101325. - pres[0][fcells[0]]) / 100.;
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf] * frac + uw_rel_perm_f[0][f] * (1 - frac);
        }
      }
    }
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) { *vo_->os() << " " << update_perm << std::endl; }
  return update_perm;
};


bool
Richards::UpdatePermeabilityDerivativeData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability derivatives?";
  if (fixed_kr_) return false;

  bool update_perm = S_->GetEvaluator(coef_key_, tag).UpdateDerivative(*S_, name_, key_, tag);
  if (update_perm) {
    const CompositeVector& drel_perm =
      S_->GetDerivative<CompositeVector>(coef_key_, tag, key_, tag);

    if (!duw_coef_key_.empty()) {
      // must also upwind
      CompositeVector& duw_rel_perm = S_->GetW<CompositeVector>(duw_coef_key_, tag, name_);
      duw_rel_perm.PutScalar(0.);

      // Upwind, only overwriting boundary faces if the wind says to do so.
      upwinding_deriv_->Update(drel_perm, duw_rel_perm, *S_);
    }
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) { *vo_->os() << " " << update_perm << std::endl; }
  return update_perm;
};


// -----------------------------------------------------------------------------
// Compute boundary condition functions at the current time.
// -----------------------------------------------------------------------------
void
Richards::ComputeBoundaryConditions_(const Tag& tag)
{
  bc_pressure_->Compute(S_->get_time(tag));
  bc_head_->Compute(S_->get_time(tag));
  bc_level_->Compute(S_->get_time(tag));
  bc_flux_->Compute(S_->get_time(tag));
  bc_seepage_->Compute(S_->get_time(tag));
  bc_seepage_infilt_->Compute(S_->get_time(tag));
}


// -----------------------------------------------------------------------------
// Push boundary conditions into the global array.
// -----------------------------------------------------------------------------
void
Richards::UpdateBoundaryConditions_(const Tag& tag, bool kr)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating BCs." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  // initialize all to 0
  for (unsigned int n = 0; n != markers.size(); ++n) {
    markers[n] = Operators::OPERATOR_BC_NONE;
    values[n] = 0.0;
  }

  // count for debugging
  std::vector<int> bc_counts;
  std::vector<std::string> bc_names;

  // Dirichlet-type boundary conditions
  // -------------------------------------
  // pressure boundary conditions -- the primary
  bc_counts.push_back(bc_pressure_->size());
  bc_names.push_back(key_);
  for (const auto& bc : *bc_pressure_) {
    int f = bc.first;
#ifdef ENABLE_DBC
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    AMANZI_ASSERT(cells.size() == 1);
#endif
    markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    values[f] = bc.second;
  }

  // head boundary conditions
  bc_counts.push_back(bc_head_->size());
  bc_names.push_back("head");
  if (bc_head_->size() > 0) {
    Errors::Message mesg("column_ID not implemented in Amanzi");
    Exceptions::amanzi_throw(mesg);

    double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);
    int z_index = mesh_->getSpaceDimension() - 1;
    const auto& gravity = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    double g = -gravity[z_index];

    for (const auto& bc : *bc_head_) {
      int f = bc.first;

      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);

      // we need to find the elevation of the surface, but finding the top edge
      // of this stack of faces is not possible currently.  The best approach
      // is instead to work with the cell.
      int col = 0;//mesh_->column_ID(cells[0]);
      double z_surf = mesh_->getFaceCentroid(mesh_->columns.getFaces(col)[0])[z_index];
      double z_wt = bc.second + z_surf;

      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      // note, here the cell centroid's z is used to relate to the column's top
      // face centroid, specifically NOT the boundary face's centroid.
      values[f] = p_atm + bc_rho_water_ * g * (z_wt - mesh_->getCellCentroid(cells[0])[z_index]);
    }
  }

  // fixed level head boundary conditions
  bc_counts.push_back(bc_head_->size());
  bc_names.push_back("head");
  if (bc_level_->size() > 0) {
    double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);
    int z_index = mesh_->getSpaceDimension() - 1;
    const auto& gravity = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    double g = -gravity[z_index];

    for (const auto& bc : *bc_level_) {
      int f = bc.first;
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = p_atm + bc_rho_water_ * g * (bc.second - mesh_->getFaceCentroid(f)[z_index]);
    }
  }

  // Neumann type boundary conditions
  // -------------------------------------
  const Epetra_MultiVector& rel_perm =
    *S_->Get<CompositeVector>(uw_coef_key_, tag).ViewComponent("face", false);

  // standard Neumann flux BCs
  bc_counts.push_back(bc_flux_->size());
  bc_names.push_back(flux_key_);

  if (!infiltrate_only_if_unfrozen_) {
    for (const auto& bc : *bc_flux_) {
      int f = bc.first;
#ifdef ENABLE_DBC
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);
#endif
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = bc.second;
      if (!kr && rel_perm[0][f] > 0.) values[f] /= rel_perm[0][f];
    }

  } else {
    // Neumann boundary conditions that turn off if temp < freezing
    const Epetra_MultiVector& temp =
      *S_->GetPtr<CompositeVector>(Keys::getKey(domain_, "temperature"), tag)
         ->ViewComponent("face");
    for (const auto& bc : *bc_flux_) {
      int f = bc.first;
#ifdef ENABLE_DBC
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);
#endif
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      if (temp[0][f] > 273.15) {
        values[f] = bc.second;
        if (!kr && rel_perm[0][f] > 0.) values[f] /= rel_perm[0][f];
      } else {
        values[f] = 0.;
      }
    }
  }

  // seepage face -- pressure <= specified value (usually 101325), outward water flux >= 0
  S_->Get<CompositeVector>(flux_key_, tag).ScatterMasterToGhosted("face");
  const Epetra_MultiVector& flux =
    *S_->Get<CompositeVector>(flux_key_, tag).ViewComponent("face", true);

  const double& p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);
  Teuchos::RCP<const CompositeVector> u = S_->GetPtr<CompositeVector>(key_, tag);
  double seepage_tol = 10.;

  bc_counts.push_back(bc_seepage_->size());
  bc_names.push_back("standard seepage");
  for (const auto& bc : *bc_seepage_) {
    int f = bc.first;
#ifdef ENABLE_DBC
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    AMANZI_ASSERT(cells.size() == 1);
#endif

    //double boundary_pressure = std::max(getFaceOnBoundaryValue(f, *u, *bc_), 101325.); // does not make sense to seep from nonsaturated cells
    double boundary_pressure =
      getFaceOnBoundaryValue(f, *u, *bc_); // does not make sense to seep from nonsaturated cells
    double boundary_flux = flux[0][f] * getBoundaryDirection(*mesh_, f);
    if (boundary_pressure > bc.second) {
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = bc.second;
    } else if (boundary_pressure < bc.second - seepage_tol) {
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = 0.;
    } else if (boundary_flux >= 0.) {
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = bc.second;
    } else {
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = 0.;
    }
  }

  // seepage face -- pressure <= p_atm, outward water flux is specified
  bc_counts.push_back(bc_seepage_infilt_->size());
  bc_names.push_back("seepage with infiltration");
  {
    Tag seepage_tag;
    if (bc_seepage_infilt_explicit_) {
      seepage_tag = tag_current_;
    } else {
      seepage_tag = tag;
    }
    const Epetra_MultiVector& flux =
      *S_->Get<CompositeVector>(flux_key_, seepage_tag).ViewComponent("face", true);
    Teuchos::RCP<const CompositeVector> u = S_->GetPtr<CompositeVector>(key_, seepage_tag);

    int i = 0;
    for (const auto& bc : *bc_seepage_infilt_) {
      int f = bc.first;
#ifdef ENABLE_DBC
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);
#endif

      double flux_seepage_tol = std::abs(bc.second) * .001;
      double boundary_pressure = getFaceOnBoundaryValue(f, *u, *bc_);
      double boundary_flux = flux[0][f] * getBoundaryDirection(*mesh_, f);
      if (i == 0)
        std::cout << "BFlux = " << boundary_flux
                  << " with constraint = " << bc.second - flux_seepage_tol << std::endl;

      if (boundary_flux < bc.second - flux_seepage_tol && boundary_pressure > p_atm + seepage_tol) {
        // both constraints are violated, either option should push things in the right direction
        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = p_atm;
        if (i == 0)
          std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux "
                    << flux[0][f] * getBoundaryDirection(*mesh_, f)
                    << " resulted in DIRICHLET pressure " << p_atm << std::endl;

      } else if (boundary_flux >= bc.second - flux_seepage_tol &&
                 boundary_pressure > p_atm - seepage_tol) {
        // max pressure condition violated
        markers[f] = Operators::OPERATOR_BC_DIRICHLET;
        values[f] = p_atm;
        if (i == 0)
          std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux "
                    << flux[0][f] * getBoundaryDirection(*mesh_, f)
                    << " resulted in DIRICHLET pressure " << p_atm << std::endl;

      } else if (boundary_flux < bc.second - flux_seepage_tol &&
                 boundary_pressure <= p_atm + seepage_tol) {
        // max infiltration violated
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = bc.second;
        if (i == 0)
          std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux "
                    << flux[0][f] * getBoundaryDirection(*mesh_, f) << " resulted in NEUMANN flux "
                    << bc.second << std::endl;

      } else if (boundary_flux >= bc.second - flux_seepage_tol &&
                 boundary_pressure <= p_atm - seepage_tol) {
        // both conditions are valid
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = bc.second;
        if (i == 0)
          std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux "
                    << flux[0][f] * getBoundaryDirection(*mesh_, f) << " resulted in NEUMANN flux "
                    << bc.second << std::endl;

      } else {
        AMANZI_ASSERT(0);
      }
      i++;
    }
  }

  // surface coupling
  bc_counts.push_back(0);
  bc_names.push_back("surface coupling (head)");

  if (coupled_to_surface_via_head_) {
    // Face is Dirichlet with value of surface head
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh("surface");
    const Epetra_MultiVector& head =
      *S_->GetPtr<CompositeVector>("surface_pressure", tag)->ViewComponent("cell", false);

    unsigned int ncells_surface = head.MyLength();
    bc_counts[bc_counts.size() - 1] = ncells_surface;

    for (unsigned int c = 0; c != ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f = surface->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
#ifdef ENABLE_DBC
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);
#endif
      // -- set that value to dirichlet
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = head[0][c];
    }
  }

  // surface coupling
  bc_counts.push_back(0);
  bc_names.push_back("surface coupling (flux)");
  if (coupled_to_surface_via_flux_) {
    // Face is Neumann with value of surface residual
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh(Keys::getDomain(ss_flux_key_));
    const Epetra_MultiVector& ss_flux =
      *S_->Get<CompositeVector>(ss_flux_key_, tag).ViewComponent("cell", false);
    unsigned int ncells_surface = ss_flux.MyLength();
    bc_counts[bc_counts.size() - 1] = ncells_surface;
    for (unsigned int c = 0; c != ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f = surface->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
#ifdef ENABLE_DBC
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);
#endif
      // -- set that value to Neumann
      markers[f] = Operators::OPERATOR_BC_NEUMANN;

      // NOTE: the flux provided by the coupler is in units of mols / s, where
      //       as Neumann BCs are in units of mols / s / A.  The right A must
      //       be chosen, as it is the subsurface mesh's face area, not the
      //       surface mesh's cell area.
      values[f] = ss_flux[0][c] / mesh_->getFaceArea(f);

      if (!kr && rel_perm[0][f] > 0.) values[f] /= rel_perm[0][f];
    }
  }

  // mark all remaining boundary conditions as zero flux conditions
  int n_default = 0;
  int nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (markers[f] == Operators::OPERATOR_BC_NONE) {
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      if (ncells == 1) {
        n_default++;
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
      }
    }
  }
  bc_names.push_back("default (zero flux)");
  bc_counts.push_back(n_default);

  // report on counts
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<int> bc_counts_global(bc_counts.size(), 0);
    mesh_->getComm()->SumAll(&bc_counts[0], &bc_counts_global[0], bc_counts.size());

    *vo_->os() << "  BCs applied:" << std::endl;

    for (int i = 0; i != bc_counts_global.size(); ++i) {
      *vo_->os() << "    " << bc_names[i] << ": " << bc_counts_global[i] << std::endl;
    }
  }
};


bool
Richards::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Modifying predictor:" << std::endl;

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);
  db_->WriteBoundaryConditions(bc_markers(), bc_values());

  // push Dirichlet data into predictor
  applyDirichletBCs(*bc_, *u->Data());

  bool changed(false);
  if (modify_predictor_bc_flux_ ||
      (modify_predictor_first_bc_flux_ && ((S_->Get<int>("cycle", Tags::DEFAULT) == 0) ||
                                           (S_->Get<int>("cycle", Tags::DEFAULT) == 1)))) {
    changed |= ModifyPredictorFluxBCs_(h, u);
  }

  if (modify_predictor_wc_) { changed |= ModifyPredictorWC_(h, u); }

  if (modify_predictor_with_consistent_faces_) { changed |= ModifyPredictorConsistentFaces_(h, u); }
  return changed;
}


bool
Richards::ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u)
{
  if (!u->Data()->HasComponent("face")) return false;

  auto& markers = bc_markers();
  auto& values = bc_values();

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

  if (flux_predictor_ == Teuchos::null) {
    flux_predictor_ =
      Teuchos::rcp(new PredictorDelegateBCFlux(S_, mesh_, matrix_diff_, wrms_, &markers, &values));
  }

  UpdatePermeabilityData_(tag_next_);
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_);

  matrix_->Init();
  matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag_next_);
  Teuchos::RCP<const CompositeVector> pres = S_->GetPtr<CompositeVector>(key_, tag_next_);
  matrix_diff_->SetDensity(rho);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
  //matrix_diff_->ApplyBCs(true, true, true);

  flux_predictor_->ModifyPredictor(h, u);
  ChangedSolution(); // mark the solution as changed, as modifying with
                     // consistent faces will then get the updated boundary
                     // conditions
  return true;
}

bool
Richards::ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications for consistent face pressures." << std::endl;

  CalculateConsistentFaces(u->Data().ptr());
  return true;
}

bool
Richards::ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u)
{
  AMANZI_ASSERT(0);
  return false;
}


// void Richards::CalculateConsistentFacesForInfiltration_(
//     const Teuchos::Ptr<CompositeVector>& u) {

//  auto& markers = bc_markers();
//  auto& values = bc_values();

//   if (vo_->os_OK(Teuchos::VERB_EXTREME))
//     *vo_->os() << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

//   if (flux_predictor_ == Teuchos::null) {
//     flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(S_next_, mesh_, matrix_,
//             wrms_, &markers, &values));
//   }

//   // update boundary conditions
//   bc_pressure_->Compute(S_->get_time(tag_next_));
//   bc_flux_->Compute(S_->get_time(tag_next_));
//   UpdateBoundaryConditions_(S_next_.ptr());

//   bool update = UpdatePermeabilityData_(S_next_.ptr());
//   Teuchos::RCP<const CompositeVector> rel_perm =
//       S_next_->GetPtr<CompositeVector>(uw_coef_key_);
//   matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
//   matrix_->CreateMFDrhsVectors();
//   Teuchos::RCP<const CompositeVector> rho = S_next_->GetPtr<CompositeVector>(mass_dens_key_);
//   Teuchos::RCP<const Epetra_Vector> gvec = S_next_->GetConstantVectorData("gravity", Tags::DEFAULT);
//   AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());
//   matrix_->ApplyBoundaryConditions(markers, values);

//   flux_predictor_->ModifyPredictor(u);
// }

void
Richards::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u)
{
  if (!u->HasComponent("face")) return; // not need

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Modifying predictor for consistent faces" << std::endl;

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

  // Using the old BCs, so should use the old rel perm?
  // update the rel perm according to the scheme of choice
  //  UpdatePermeabilityData_(tag_next_);
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_);

  S_->GetEvaluator(mass_dens_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag_next_);

  // Update the preconditioner with darcy and gravity fluxes
  matrix_->Init();
  matrix_diff_->SetDensity(rho);
  matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, u);
  matrix_diff_->ApplyBCs(true, true, true);

  // derive the consistent faces, involves a solve
  db_->WriteVector(" p_cf guess:", u.ptr(), true);
  matrix_diff_->UpdateConsistentFaces(*u);
  db_->WriteVector(" p_cf soln:", u.ptr(), true);
}


// /* ******************************************************************
// * Clip pressure using pressure threshold.
// ****************************************************************** */
// void Richards::ClipHydrostaticPressure(double pmin, Epetra_MultiVector& p)
// {
//   int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
//   for (int c = 0; c < ncells_owned; c++) p[0][c] = std::max(p[0][c], pmin);
// }


// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool
Richards::IsAdmissible(Teuchos::RCP<const TreeVector> up)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Checking admissibility..." << std::endl;

  // For some reason, wandering PKs break most frequently with an unreasonable
  // pressure.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> pres = up->Data();
  double minT, maxT;

  const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
  double minT_c(1.e15), maxT_c(-1.e15);
  int min_c(-1), max_c(-1);
  for (int c = 0; c != pres_c.MyLength(); ++c) {
    if (pres_c[0][c] < minT_c) {
      minT_c = pres_c[0][c];
      min_c = c;
    }
    if (pres_c[0][c] > maxT_c) {
      maxT_c = pres_c[0][c];
      max_c = c;
    }
  }

  double minT_f(1.e15), maxT_f(-1.e15);
  int min_f(-1), max_f(-1);
  if (pres->HasComponent("face")) {
    const Epetra_MultiVector& pres_f = *pres->ViewComponent("face", false);
    for (int f = 0; f != pres_f.MyLength(); ++f) {
      if (pres_f[0][f] < minT_f) {
        minT_f = pres_f[0][f];
        min_f = f;
      }
      if (pres_f[0][f] > maxT_f) {
        maxT_f = pres_f[0][f];
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
    *vo_->os() << "    Admissible p? (min/max): " << minT << ",  " << maxT << std::endl;
  }

  if (minT < -1.e9 || maxT > 1.e8) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << " is not admissible, as it is not within bounds of constitutive models:"
                 << std::endl;

      Teuchos::RCP<const Comm_type> comm_p = mesh_->getComm();
      Teuchos::RCP<const MpiComm_type> mpi_comm_p =
        Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
      const MPI_Comm& comm = mpi_comm_p->Comm();

      ENorm_t global_minT_c, local_minT_c;
      ENorm_t global_maxT_c, local_maxT_c;

      local_minT_c.value = minT_c;
      local_minT_c.gid = pres_c.Map().GID(min_c);
      local_maxT_c.value = maxT_c;
      local_maxT_c.gid = pres_c.Map().GID(max_c);

      MPI_Allreduce(&local_minT_c, &global_minT_c, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
      MPI_Allreduce(&local_maxT_c, &global_maxT_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
      *vo_->os() << "   cells (min/max): [" << global_minT_c.gid << "] " << global_minT_c.value
                 << ", [" << global_maxT_c.gid << "] " << global_maxT_c.value << std::endl;

      if (pres->HasComponent("face")) {
        const Epetra_MultiVector& pres_f = *pres->ViewComponent("face", false);
        ENorm_t global_minT_f, local_minT_f;
        ENorm_t global_maxT_f, local_maxT_f;

        local_minT_f.value = minT_f;
        local_minT_f.gid = pres_f.Map().GID(min_f);
        local_maxT_f.value = maxT_f;
        local_maxT_f.gid = pres_f.Map().GID(max_f);

        MPI_Allreduce(&local_minT_f, &global_minT_f, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
        MPI_Allreduce(&local_maxT_f, &global_maxT_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
        *vo_->os() << "   faces (min/max): [" << global_minT_f.gid << "] " << global_minT_f.value
                   << ", [" << global_maxT_f.gid << "] " << global_maxT_f.value << std::endl;
      }
    }
    return false;
  }
  return true;
}


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Richards::ModifyCorrection(double h,
                           Teuchos::RCP<const TreeVector> res,
                           Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // if the primary variable has boundary face, this is for upwinding rel
  // perms and is never actually used.  Make sure it does not go to undefined
  // pressures.
  if (du->Data()->HasComponent("boundary_face")) {
    du->Data()->ViewComponent("boundary_face")->PutScalar(0.);
  }

  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp = du->Data()->begin(); comp != du->Data()->end();
       ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp, false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2
                 << std::endl;
    }
  }

  // limit by capping corrections when they cross atmospheric pressure
  // (where pressure derivatives are discontinuous)
  int my_limited = 0;
  int n_limited_spurt = 0;
  if (patm_limit_ > 0.) {
    double patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);
    for (CompositeVector::name_iterator comp = du->Data()->begin(); comp != du->Data()->end();
         ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp, false);
      const Epetra_MultiVector& u_c = *u->Data()->ViewComponent(*comp, false);

      for (int c = 0; c != du_c.MyLength(); ++c) {
        if ((u_c[0][c] < patm) && (u_c[0][c] - du_c[0][c] > patm + patm_limit_)) {
          du_c[0][c] = u_c[0][c] - (patm + patm_limit_);
          my_limited++;
        } else if ((u_c[0][c] > patm) && (u_c[0][c] - du_c[0][c] < patm - patm_limit_)) {
          du_c[0][c] = u_c[0][c] - (patm - patm_limit_);
          my_limited++;
        }
      }
    }
    mesh_->getComm()->MaxAll(&my_limited, &n_limited_spurt, 1);
  }

  if (n_limited_spurt > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) { *vo_->os() << "  limiting the spurt." << std::endl; }
  }

  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp = du->Data()->begin(); comp != du->Data()->end();
       ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp, false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2
                 << std::endl;
    }
  }

  // Limit based on a max pressure change
  my_limited = 0;
  int n_limited_change = 0;
  if (p_limit_ >= 0.) {
    for (CompositeVector::name_iterator comp = du->Data()->begin(); comp != du->Data()->end();
         ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp, false);

      double max;
      du_c.NormInf(&max);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Max pressure correction (" << *comp << ") = " << max << std::endl;
      }

      for (int c = 0; c != du_c.MyLength(); ++c) {
        if (std::abs(du_c[0][c]) > p_limit_) {
          du_c[0][c] = ((du_c[0][c] > 0) - (du_c[0][c] < 0)) * p_limit_;
          my_limited++;
        }
      }
    }

    mesh_->getComm()->MaxAll(&my_limited, &n_limited_change, 1);
  }

  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp = du->Data()->begin(); comp != du->Data()->end();
       ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp, false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2
                 << std::endl;
    }
  }

  if (n_limited_change > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) { *vo_->os() << "  limited by pressure." << std::endl; }
  }

  if (n_limited_spurt > 0) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING;
  } else if (n_limited_change > 0) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  }

  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

} // namespace Flow
} // namespace Amanzi
