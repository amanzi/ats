/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
#include "Reductions.hh"
#include "Point.hh"
#include "DataStructuresHelpers.hh"

#include "EvaluatorModelCVByMaterial.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionWithGravity.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Accumulation.hh"

#include "PK_Helpers.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
//#include "upwind_gravity_flux.hh"

#include "predictor_delegate_bc_flux.hh"
//#include "BoundaryFlux.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

const std::string Richards::pk_type_ = "richards flow";


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
Richards::Richards(const Comm_ptr_type& comm,
                   Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S)
  : PK_PhysicalBDF_Default(comm, pk_tree, glist, S),
    coupled_to_surface_via_head_(false),
    coupled_to_surface_via_flux_(false),
    infiltrate_only_if_unfrozen_(false),
    modify_predictor_with_consistent_faces_(false),
    modify_predictor_wc_(false),
    modify_predictor_bc_flux_(false),
    modify_predictor_first_bc_flux_(false),
    upwind_from_prev_flux_(false),
    clobber_boundary_flux_dir_(false),
    // perm_scale_(1.),
    jacobian_(false),
    jacobian_lag_(0),
    iter_(0),
    iter_counter_time_(0.),
    perm_scale_(1.),
    fixed_kr_(false)
{}


void
Richards::modifyParameterList()
{
  // set some defaults for inherited PKs
  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set<std::string>("conserved quantity key suffix", "water_content");

  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .5 * .1 * 55000.); // phi * s * nl

  // scaling for permeability for better "nondimensionalization"
  double perm_scale = plist_->get<double>("permeability rescaling", 1.e7);
  S_->ConstantsList().sublist("permeability_rescaling").set<double>("value", perm_scale);

  PK_PhysicalBDF_Default::modifyParameterList();
}


void
Richards::parseParameterList()
{
  // parse inherited lists
  PK_PhysicalBDF_Default::parseParameterList();

  // get field names
  mass_dens_key_ = Keys::readKey(*plist_, domain_, "mass density", "mass_density_liquid");
  molar_dens_key_ = Keys::readKey(*plist_, domain_, "molar density", "molar_density_liquid");
  perm_key_ = Keys::readKey(*plist_, domain_, "permeability", "permeability");
  coef_key_ = Keys::readKey(
    *plist_, domain_, "relative hydraulic conductivity", "relative_hydraulic_conductivity");
  uw_coef_key_ = Keys::readKey(
    *plist_, domain_, "upwinded conductivity", "upwind_relative_hydraulic_conductivity");
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

  // source terms
  is_source_term_ = plist_->get<bool>("source term", false);
  if (is_source_term_) {
    if (source_key_.empty()) {
      source_key_ = Keys::readKey(*plist_, domain_, "source", "water_source");
    }
    source_term_is_differentiable_ = plist_->get<bool>("source term is differentiable", true);
    explicit_source_ = plist_->get<bool>("explicit source term", false);
  }

  // coupling to surface
  coupled_to_surface_via_flux_ = plist_->get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    Key domain_surf = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");
    ss_flux_key_ =
      Keys::readKey(*plist_, domain_surf, "surface-subsurface flux", "surface_subsurface_flux");
  }

  coupled_to_surface_via_head_ = plist_->get<bool>("coupled to surface via head", false);
  if (coupled_to_surface_via_head_) {
    Key domain_surf = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");
    ss_primary_key_ = Keys::readKey(*plist_, domain_surf, "pressure", "pressure");
  }

  // boundary conditions
  auto bc_plist = Teuchos::sublist(plist_, "boundary conditions");
  std::vector<std::string> bcs_found;

  for (auto& bc_sublist : *bc_plist) {
    if (bc_sublist.first == "pressure") {
      // -- Dirichlet conditions
      std::string bc_pressure_name = name_ + "_bcs_pressure";
      Teuchos::ParameterList& bc_pressure_list = S_->GetEvaluatorList(bc_pressure_name);
      bc_pressure_list.set("pressure", bc_plist->sublist("pressure"));
      bc_pressure_list.set<std::string>("evaluator type", "independent variable patch")
        .set<std::string>("function list name", "pressure")
        .set<std::string>("function inner list name", "boundary pressure [Pa]");
      bcs_found.emplace_back(bc_pressure_name);

    } else if (bc_sublist.first == "water flux") {
      // -- Neumann conditions
      std::string bc_flux_name = name_ + "_bcs_water_flux";
      Teuchos::ParameterList& bc_flux_list = S_->GetEvaluatorList(bc_flux_name);
      bc_flux_list.set("water flux", bc_plist->sublist("water flux"));
      bc_flux_list.set<std::string>("evaluator type", "independent variable patch")
        .set<std::string>("function list name", "water flux")
        .set<std::string>("function inner list name", "outward water flux [mol m^-2 s^-1]");
      bcs_found.emplace_back(bc_flux_name);
      // } else if (bc_sublist.first == "seepage face pressure") {
      // requires the addition of an EvaluatorFlowBCsSeepage
      //   // -- seepage condition where a max pressure is supplied
      //   std::string bc_pressure_name = name_ + "_bcs_seepage_face_pressure";
      //   Teuchos::ParameterList& bc_pressure_list = S_->GetEvaluatorList(bc_pressure_name);
      //   bc_pressure_list.set("pressure", bc_plist->sublist("pressure"));
      //   bc_pressure_list.set<std::string>("evaluator type", "independent variable patch")
      //     .set<std::string>("function list name", "pressure")
      //     .set<std::string>("function inner list name", "boundary pressure [Pa]");
      //   bcs_found.emplace_back(bc_pressure_name);
    }
  }

  if (coupled_to_surface_via_flux_) {
    // Neumann condition, given by the vector
    std::string bc_flux_name = name_ + "_bcs_surface_coupling_via_flux";
    Teuchos::ParameterList& bc_flux_list = S_->GetEvaluatorList(bc_flux_name);
    bc_flux_list.set<std::string>("evaluator type", "vector as patch")
      .set<Teuchos::Array<std::string>>("dependencies", std::vector<std::string>{ ss_flux_key_ })
      .set<std::string>("region", "surface");
    bcs_found.emplace_back(bc_flux_name);
  }

  // and the boundary condition aggregator
  S_->GetEvaluatorList(name_ + "_bcs")
    .set<std::string>("evaluator type", "boundary condition aggregator")
    .set<std::string>("entity kind", "face")
    .set<Teuchos::Array<std::string>>("dependencies", bcs_found);
}

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void
Richards::setup()
{
  PK_PhysicalBDF_Default::setup();
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
  // is dynamic mesh?  If so, get a key for indicating when the mesh has changed.
  if (!deform_key_.empty()) S_->RequireEvaluator(deform_key_, tag_next_);

  S_->Require<double>("permeability_rescaling", Tags::DEFAULT);

  //
  // Diffusion Operators
  // ------------------------------------------------------------------
  // -- Boundary conditions.
  auto& bc_fac = S_->Require<Operators::BCs, Operators::BCs_Factory>(name_ + "_bcs", tag_next_);
  bc_fac.set_mesh(mesh_);
  bc_fac.set_entity_kind(AmanziMesh::Entity_kind::FACE);
  bc_fac.set_dof_type(WhetStone::DOF_Type::SCALAR);
  S_->RequireEvaluator(name_ + "_bcs", tag_next_);

  if (plist_->sublist("boundary conditions").isSublist("pressure"))
    S_->Require<MultiPatch<double>, MultiPatchSpace>(name_ + "_bcs_pressure", tag_next_)
      .set_flag(Operators::OPERATOR_BC_DIRICHLET);
  if (plist_->sublist("boundary conditions").isSublist("water flux"))
    S_->Require<MultiPatch<double>, MultiPatchSpace>(name_ + "_bcs_water_flux", tag_next_)
      .set_flag(Operators::OPERATOR_BC_NEUMANN);
  if (plist_->sublist("boundary conditions").isSublist("seepage face pressure"))
    S_->Require<MultiPatch<double>, MultiPatchSpace>(name_ + "_bcs_seepage_face_pressure",
                                                     tag_next_)
      .set_flag(Operators::OPERATOR_BC_CONDITIONAL);
  if (coupled_to_surface_via_flux_) {
    auto& mps = S_->Require<MultiPatch<double>, MultiPatchSpace>(
      name_ + "_bcs_surface_coupling_via_flux", tag_next_);
    mps.set_flag(Operators::OPERATOR_BC_NEUMANN);
    // THIS IS BROKEN AND FRAGILE.... HOW DO WE DO THIS CORRECTLY? --ETC
    mps.addPatch("surface", AmanziMesh::Entity_kind::FACE, 1);
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
  AMANZI_ASSERT(clobber_boundary_flux_dir_ == false);

  // is dynamic mesh?  If so, get a key for indicating when the mesh has changed.
  if (!deform_key_.empty()) S_->RequireEvaluator(deform_key_, tag_next_);

  // what upwinding method to use
  std::string method_name =
    plist_->get<std::string>("relative permeability method", "upwind with Darcy flux");
  // if (method_name == "upwind with gravity") {
  //   upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_, tag_next_, K_));
  //   Krel_method_ = Operators::UPWIND_METHOD_GRAVITY;
  // } else if (method_name == "cell centered") {
  if (method_name == "cell centered") {
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
  PKHelpers::requireNonlinearDiffusionCoefficient(uw_coef_key_, tag_next_, coef_location, *S_);

  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  mfd_plist.set("gravity", true);

  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.CreateWithGravity(mfd_plist, mesh_);
  matrix_ = matrix_diff_->global_operator();

  // -- create the operator, data for flux directions
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = opfactory.CreateWithGravity(face_diff_list, mesh_);

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
  preconditioner_diff_ = opfactory.CreateWithGravity(mfd_pc_plist, mesh_);
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
      PKHelpers::requireNonlinearDiffusionCoefficient(
        duw_coef_key_, tag_next_, "upwind: face", *S_);

      // note, this is here to be consistent -- unclear whether the 1.e-3 is useful or not?
      double flux_eps = plist_->get<double>("upwind flux epsilon", 1.e-5);
      upwinding_deriv_ = Teuchos::rcp(
        new Operators::UpwindTotalFlux(name_, tag_next_, flux_dir_key_, 1.e-3 * flux_eps));
    } else {
      // FV -- no upwinding of derivative
      duw_coef_key_ = std::string();
    }
  }

  //
  // Accumluation Operator
  // ------------------------------------------------------------------
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set<std::string>("entity kind", "cell");
  preconditioner_acc_ =
    Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  //
  // Source term
  // ------------------------------------------------------------------
  if (is_source_term_) {
    PKHelpers::requireAtNext(source_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    if (source_term_is_differentiable_) {
      // require derivative of source
      S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        source_key_, tag_next_, key_, tag_next_);
    }
  }

  //
  // Coupling to surface
  // ------------------------------------------------------------------
  // -- coupling done by a Neumann condition
  if (coupled_to_surface_via_flux_) {
    S_->Require<CompositeVector, CompositeVectorSpace>(ss_flux_key_, tag_next_)
      .SetMesh(S_->GetMesh(Keys::getDomain(ss_flux_key_)))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- coupling done by a Dirichlet condition
  if (coupled_to_surface_via_head_) {
    S_->Require<CompositeVector, CompositeVectorSpace>(ss_primary_key_, tag_next_)
      .SetMesh(S_->GetMesh(Keys::getDomain(ss_primary_key_)))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- Make sure coupling isn't flagged multiple ways.
  if (coupled_to_surface_via_flux_ && coupled_to_surface_via_head_) {
    Errors::Message message("Richards PK requested both flux and head coupling -- choose one.");
    Exceptions::amanzi_throw(message);
  }

  //
  // Require fields and evaluators
  // ------------------------------------------------------------------

  // -- pressure, the primary variable
  //  NOTE: no need to require evaluator for p here, either at the old or new
  //  times, as this was done in pk_physical.  All we have to do is set the
  //  structure.
  compute_boundary_values_ = plist_->get<bool>("compute boundary values", false);
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_)
    .Update(*matrix_->getRangeMap())
    ->SetGhosted()
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

  // -- flux is managed here as a primary variable
  PKHelpers::requireAtNext(flux_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // -- also need a velocity, but only for vis/diagnostics, so might as well
  // -- only keep at NEXT
  PKHelpers::requireAtNext(velocity_key_, Tags::NEXT, *S_, name_)
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

  // NOT YET IMPLEMENTED!
  AMANZI_ASSERT(!modify_predictor_with_consistent_faces_);
  AMANZI_ASSERT(!modify_predictor_first_bc_flux_);
  AMANZI_ASSERT(!modify_predictor_wc_);

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
  // -- permeability tensor
  auto& tv_fac = S_->Require<TensorVector, TensorVector_Factory>(perm_key_, tag_next_);
  tv_fac.getMap().SetMesh(mesh_)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(perm_key_, tag_next_);

  // -- water content, and evaluator, and derivative for PC
  PKHelpers::requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    conserved_key_, tag_next_, key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  PKHelpers::requireAtCurrent(conserved_key_, tag_current_, *S_, name_);

  // -- Water retention evaluators
  // -- saturation
  PKHelpers::requireAtNext(sat_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  PKHelpers::requireAtNext(sat_gas_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  auto& wrm = S_->RequireEvaluator(sat_key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  PKHelpers::requireAtCurrent(sat_key_, tag_current_, *S_, name_);

  // -- rel perm
  PKHelpers::requireAtNext(coef_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

  // -- molar density used to infer liquid Darcy velocity from flux
  PKHelpers::requireAtNext(molar_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // -- liquid mass density for the gravity fluxes
  PKHelpers::requireAtNext(mass_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void
Richards::initialize()
{
  // Initialize via hydrostatic balance
  if (!S_->GetRecordW(key_, tag_next_, name_).initialized()) InitializeHydrostatic_(tag_next_);

  // Initialize in the standard ways
  PK_PhysicalBDF_Default::initialize();

  // Set extra fields as initialized -- these don't currently have evaluators,
  // and will be initialized in the call to commit_state()
  S_->GetW<CompositeVector>(uw_coef_key_, tag_next_, uw_coef_key_).putScalar(1.0);
  S_->GetRecordW(uw_coef_key_, tag_next_, uw_coef_key_).set_initialized();

  if (!duw_coef_key_.empty()) {
    S_->GetW<CompositeVector>(duw_coef_key_, tag_next_, duw_coef_key_).putScalar(1.0);
    S_->GetRecordW(duw_coef_key_, tag_next_, duw_coef_key_).set_initialized();
  }

  S_->GetW<CompositeVector>(flux_key_, tag_next_, getName()).putScalar(0.0);
  S_->GetRecordW(flux_key_, tag_next_, getName()).set_initialized();
  PKHelpers::changedEvaluatorPrimary(flux_key_, tag_next_, *S_);

  S_->GetW<CompositeVector>(flux_dir_key_, tag_next_, getName()).putScalar(0.0);
  S_->GetRecordW(flux_dir_key_, tag_next_, getName()).set_initialized();
  S_->GetW<CompositeVector>(velocity_key_, Tags::NEXT, getName()).putScalar(0.0);
  S_->GetRecordW(velocity_key_, Tags::NEXT, getName()).set_initialized();

  // operators
  const auto& g = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  const auto& bc = S_->GetPtr<Operators::BCs>(name_ + "_bcs", tag_next_);
  matrix_diff_->SetGravity(g);
  matrix_diff_->SetBCs(bc, bc);

  S_->GetEvaluator(perm_key_, Tags::DEFAULT).Update(*S_, name_);
  auto K = S_->GetPtr<TensorVector>(perm_key_, Tags::DEFAULT);
  matrix_diff_->SetTensorCoefficient(K);

  preconditioner_diff_->SetGravity(g);
  preconditioner_diff_->SetBCs(bc, bc);
  preconditioner_diff_->SetTensorCoefficient(K);

  face_matrix_diff_->SetGravity(g);
  face_matrix_diff_->SetBCs(bc, bc);
  face_matrix_diff_->SetTensorCoefficient(K);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
};


void
Richards::InitializeHydrostatic_(const Tag& tag)
{
  // constant head over the surface
  if (plist_->sublist("initial condition").isParameter("hydrostatic head [m]")) {
    double head_wt = plist_->sublist("initial condition").get<double>("hydrostatic head [m]");
    double rho =
      plist_->sublist("initial condition").get<double>("hydrostatic water density [kg m^-3]");

    const auto& gvec = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    int z_index = mesh_->getSpaceDimension() - 1;
    double g = -gvec[z_index];
    double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    // set pressure on the column of faces and cells
    Teuchos::RCP<CompositeVector> pres = S_->GetPtrW<CompositeVector>(key_, tag, getName());
    bool has_faces = pres->hasComponent("face");
    Kokkos::View<int*> touched("touched", 0);

    { // context for viewcomponent -- do cells
      auto pres_c = pres->viewComponent("cell", false);
      decltype(pres_c) pres_f;
      if (has_faces) {
        pres_f = pres->viewComponent("face", false);
        Kokkos::resize(touched, pres_f.extent(0));
        Kokkos::deep_copy(touched, 0);
      }

      AMANZI_ASSERT(mesh_->columns->num_columns_owned >= 0);
      const AmanziMesh::MeshCache& m = mesh_->getCache();

      Kokkos::parallel_for(
        "Richards::InitializeHydrostatic cells",
        m.columns.num_columns_owned,
        KOKKOS_LAMBDA(const int col) {
          const auto& col_cells = m.columns.getCells<MemSpace_kind::DEVICE>(col);
          const auto& col_faces = m.columns.getFaces<MemSpace_kind::DEVICE>(col);
          double z_wt = m.getFaceCentroid(col_faces(0))[z_index] + head_wt;

          if (has_faces) {
            pres_f(col_faces(0), 0) = p_atm + rho * g * head_wt;
            touched(col_faces(0)) = 1;
          }

          for (int lcv_c = 0; lcv_c != col_cells.size(); ++lcv_c) {
            AmanziMesh::Entity_ID c = col_cells(lcv_c);
            AmanziMesh::Entity_ID f = col_faces(lcv_c + 1);
            pres_c(c, 0) = p_atm + rho * g * (z_wt - m.getCellCentroid(c)[z_index]);
            if (has_faces) {
              pres_f(f, 0) = p_atm + rho * g * (z_wt - m.getFaceCentroid(f)[z_index]);
              touched(f) = 1;
            }
          }
        });
    }

    if (has_faces) {
      // communicate, then deal with horizontal-normal faces
      pres->scatterMasterToGhosted("cell");
      {
        auto pres_c = pres->viewComponent("cell", false);
        auto pres_f = pres->viewComponent("face", false);
        const AmanziMesh::MeshCache& m = mesh_->getCache();
        Kokkos::parallel_for(
          "Richards::InitializeHydrostatic faces", pres_f.extent(0), KOKKOS_LAMBDA(const int f) {
            if (!touched(f)) {
              auto f_cells = m.getFaceCells(f);
              if (f_cells.size() == 1) {
                // boundary face, use the cell value as the water table is
                // assumed to parallel the cell structure
                pres_f(f, 0) = pres_c(f_cells(0), 0);
                touched(f) = 1;
              } else {
                // interpolate between cells
                pres_f(f, 0) = (pres_c(f_cells(0), 0) + pres_c(f_cells(1), 0)) / 2.;
                touched(f) = 1;
              }
            }
          });
      }
    }
    S_->GetRecordW(key_, tag, getName()).set_initialized();
  }

  // constant head datum
  if (plist_->sublist("initial condition").isParameter("hydrostatic water level [m]")) {
    AMANZI_ASSERT(false);
    //     double z_wt = plist_->sublist("initial condition").get<double>("hydrostatic water level [m]");
    //     double rho =
    //       plist_->sublist("initial condition").get<double>("hydrostatic water density [kg m^-3]");

    //     int z_index = mesh_->getSpaceDimension() - 1;
    //     const auto& gravity = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    //     double g = -gravity[z_index];

    //     double p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    //     Teuchos::RCP<CompositeVector> pres = S_->GetPtrW<CompositeVector>(key_, tag, getName());
    //     Epetra_MultiVector& pres_c = *pres->viewComponent("cell", false);
    //     for (int c = 0; c != pres_c.MyLength(); ++c) {
    //       pres_c[0][c] = p_atm + rho * g * (z_wt - mesh_->getCellCentroid(c)[z_index]);
    //     }

    //     if (pres->hasComponent("face")) {
    //       Epetra_MultiVector& pres_f = *pres->viewComponent("face", false);
    //       for (int f = 0; f != pres_f.MyLength(); ++f) {
    //         pres_f[0][f] = p_atm + rho * g * (z_wt - mesh_->getFaceCentroid(f)[z_index]);
    //       }
    //     }
    //     S_->GetRecordW(key_, tag, getName()).set_initialized();
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
Richards::commitStep(double t_old, double t_new, const Tag& tag_next)
{
  // saves primary variable, conserved quantity
  PK_PhysicalBDF_Default::commitStep(t_old, t_new, tag_next);

  // also save saturation
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(sat_key_, tag_current, tag_next, *S_);
  if (S_->HasRecordSet(sat_ice_key_)) {
    PKHelpers::assign(sat_ice_key_, tag_current, tag_next, *S_);
  }
};


// -----------------------------------------------------------------------------
// Check for controls on saturation
// -----------------------------------------------------------------------------
bool
Richards::isValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Validating time step." << std::endl;

  if (sat_change_limit_ > 0.0) {
    const auto& sl_new = *S_->Get<CompositeVector>(sat_key_, tag_next_).getComponent("cell", false);
    const auto& sl_old =
      *S_->Get<CompositeVector>(sat_key_, tag_current_).getComponent("cell", false);
    MultiVector_type dsl(sl_new.getMap(), 1);
    dsl.update(-1., sl_old, 1., sl_new, 0.);
    auto change = Reductions::reduceAllMaxLoc(*dsl.getVector(0));

    if (change.val > sat_change_limit_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max sl change=" << change.val
                   << " > limit=" << sat_change_limit_ << " at cell GID " << change.loc
                   << std::endl;
      return false;
    }
  }

  if (S_->HasRecordSet(sat_ice_key_) && (sat_ice_change_limit_ > 0.0)) {
    const Vector_type& si_new = *S_->GetPtr<CompositeVector>(sat_ice_key_, tag_next_)
                                   ->getComponent("cell", false)
                                   ->getVector(0);
    const Vector_type& si_old = *S_->GetPtr<CompositeVector>(sat_ice_key_, tag_current_)
                                   ->getComponent("cell", false)
                                   ->getVector(0);
    Vector_type dsi(si_new);
    dsi.update(-1., si_old, 1.);
    auto change = Reductions::reduceAllMaxLoc(dsi);

    if (change.val > sat_ice_change_limit_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max si change=" << change.val
                   << " > limit=" << sat_ice_change_limit_ << " at cell GID " << change.loc
                   << std::endl;
      return false;
    }
  }
  return PK_PhysicalBDF_Default::isValidStep();
}


// -----------------------------------------------------------------------------
// Update any diagnostic variables prior to vis (in this case velocity field).
// -----------------------------------------------------------------------------
void
Richards::calculateDiagnostics(const Tag& tag)
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

  // derive fluxes
  Teuchos::RCP<CompositeVector> flux = S_->GetPtrW<CompositeVector>(flux_key_, tag_next_, name_);
  matrix_diff_->UpdateFlux(pres.ptr(), flux.ptr());
  //  UpdateVelocity_(tag);
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

  bool update_perm = S_->GetEvaluator(coef_key_, tag).Update(*S_, name_);
  Teuchos::RCP<const CompositeVector> rel_perm = S_->GetPtr<CompositeVector>(coef_key_, tag);

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
          S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " flux dir")) {
        S_->GetEvaluator(perm_key_, Tags::DEFAULT).Update(*S_, name_ + " flux_dir");
        auto K = S_->GetPtr<TensorVector>(perm_key_, Tags::DEFAULT);
        face_matrix_diff_->SetTensorCoefficient(K);
      }
      face_matrix_diff_->SetDensity(rho);
      face_matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
      face_matrix_diff_->UpdateFlux(pres.ptr(), flux_dir.ptr());

      // if (clobber_boundary_flux_dir_) {
      //   Epetra_MultiVector& flux_dir_f = *flux_dir->viewComponent("face", false);

      //   auto& markers = bc_markers();
      //   auto& values = bc_values();

      //   for (int f = 0; f != markers.size(); ++f) {
      //     if (markers[f] == Operators::OPERATOR_BC_NEUMANN) {
      //       AmanziMesh::Entity_ID_List cells;
      //       cells = mesh_->getFaceCells(f);
      //       AMANZI_ASSERT(cells.size() == 1);
      //       int c = cells[0];
      //       AmanziMesh::Entity_ID_List faces;
      //       std::vector<int> dirs;
      //       mesh_->getCellFacesAndDirections(c, &faces, &dirs);
      //       int i = std::find(faces.begin(), faces.end(), f) - faces.begin();

      //       flux_dir_f[0][f] = values[f] * dirs[i];
      //     }
      //   }
      // }
    }

    update_perm |= update_dir;
  }

  if (update_perm) {
    Teuchos::RCP<CompositeVector> uw_rel_perm =
      S_->GetPtrW<CompositeVector>(uw_coef_key_, tag, uw_coef_key_);

    // Move rel perm on boundary_faces into uw_rel_perm on faces
    const auto& vandelay = mesh_->getBoundaryFaceImporter();
    const auto& rel_perm_bf = *rel_perm->getComponent("boundary_face", false);
    {
      auto& uw_rel_perm_f = *uw_rel_perm->getComponent("face", false);
      uw_rel_perm_f.doExport(rel_perm_bf, vandelay, Tpetra::CombineMode::INSERT);
    }

    // Upwind, only overwriting boundary faces if the wind says to do so.
    upwinding_->Update(*rel_perm, *uw_rel_perm, *S_);

    if (clobber_policy_ == "clobber") {
      auto& uw_rel_perm_f = *uw_rel_perm->getComponent("face", false);
      uw_rel_perm_f.doExport(rel_perm_bf, vandelay, Tpetra::CombineMode::INSERT);
    } else if (clobber_policy_ == "max") {
      AMANZI_ASSERT(false);
      // Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->viewComponent("face", false);
      // const auto& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
      // const auto& bfmap = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
      // for (int bf = 0; bf != rel_perm_bf.MyLength(); ++bf) {
      //   auto f = fmap.LID(bfmap.GID(bf));
      //   if (rel_perm_bf[0][bf] > uw_rel_perm_f[0][f]) { uw_rel_perm_f[0][f] = rel_perm_bf[0][bf]; }
      // }
    } else if (clobber_policy_ == "unsaturated") {
      // clobber only when the interior cell is unsaturated
      AMANZI_ASSERT(false);
      // Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->viewComponent("face", false);
      // const Epetra_MultiVector& pres =
      //   *S_->Get<CompositeVector>(key_, tag).viewComponent("cell", false);
      // const auto& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
      // const auto& bfmap = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
      // for (int bf = 0; bf != rel_perm_bf.MyLength(); ++bf) {
      //   auto f = fmap.LID(bfmap.GID(bf));
      //   AmanziMesh::Entity_ID_List fcells;
      //   fcells = mesh_->getFaceCells(f);
      //   AMANZI_ASSERT(fcells.size() == 1);
      //   if (pres[0][fcells[0]] < 101225.) {
      //     uw_rel_perm_f[0][f] = rel_perm_bf[0][bf];
      //   } else if (pres[0][fcells[0]] < 101325.) {
      //     double frac = (101325. - pres[0][fcells[0]]) / 100.;
      //     uw_rel_perm_f[0][f] = rel_perm_bf[0][bf] * frac + uw_rel_perm_f[0][f] * (1 - frac);
      //   }
      // }
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
      CompositeVector& duw_rel_perm = S_->GetW<CompositeVector>(duw_coef_key_, tag, duw_coef_key_);
      duw_rel_perm.putScalar(0.);

      // Upwind, only overwriting boundary faces if the wind says to do so.
      upwinding_deriv_->Update(drel_perm, duw_rel_perm, *S_);
    }
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) { *vo_->os() << " " << update_perm << std::endl; }
  return update_perm;
};


// -----------------------------------------------------------------------------
// Push boundary conditions into the global array.
// -----------------------------------------------------------------------------
void
Richards::UpdateBoundaryConditions_(const Tag& tag, bool kr)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating BCs." << std::endl;

  // master does this in the course of normal UpdateBoundaryConditions (for
  // seepage faces) -- do we rely on it?
  S_->Get<CompositeVector>(flux_key_, tag).scatterMasterToGhosted("face");

  S_->GetEvaluator(name_ + "_bcs", tag_next_).Update(*S_, name_);
};


bool
Richards::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Modifying predictor:" << std::endl;

  // update boundary conditions
  UpdateBoundaryConditions_(tag_next_);
  // db_->WriteBoundaryConditions(bc_->model(), bc_->value());

  // push Dirichlet data into predictor
  const auto& bcs = S_->Get<Operators::BCs>(name_ + "_bcs", tag_next_);
  PKHelpers::applyDirichletBCs(bcs, *u->getData());

  bool changed(false);
  if (modify_predictor_bc_flux_ ||
      (modify_predictor_first_bc_flux_ && ((S_->Get<int>("cycle", Tags::DEFAULT) == 0) ||
                                           (S_->Get<int>("cycle", Tags::DEFAULT) == 1)))) {
    changed |= ModifyPredictorFluxBCs_(h, u);
  }

  // if (modify_predictor_wc_) { changed |= ModifyPredictorWC_(h, u); }

  // if (modify_predictor_with_consistent_faces_) { changed |= ModifyPredictorConsistentFaces_(h, u); }
  return changed;
}


bool
Richards::ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u)
{
  if (!u->getData()->hasComponent("face")) return false;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

  if (flux_predictor_ == Teuchos::null) {
    const auto& bcs = S_->Get<Operators::BCs>(name_ + "_bcs", tag_next_);

    auto wrm_eval_as_eval = S_->GetEvaluatorPtr(sat_key_, tag_next_);
    auto wrm_eval_as_wrm =
      Teuchos::rcp_dynamic_cast<PredictorDelegateBCFlux::WRMEval_type>(wrm_eval_as_eval);
    if (wrm_eval_as_wrm == Teuchos::null) {
      Errors::Message msg("To use Richards option \"modify predictor for flux BCs\", may only use "
                          "WRM evaluator of type \"wrm van Genuchten\"");
      Exceptions::amanzi_throw(msg);
    }
    flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(
      S_, mesh_, matrix_diff_, wrm_eval_as_wrm->getModels(), bcs.model(), bcs.value()));
  }

  UpdatePermeabilityData_(tag_next_);
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_);

  matrix_->Zero();
  matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag_next_);
  Teuchos::RCP<const CompositeVector> pres = S_->GetPtr<CompositeVector>(key_, tag_next_);
  matrix_diff_->SetDensity(rho);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
  //matrix_diff_->ApplyBCs(true, true, true);

  flux_predictor_->ModifyPredictor(h, u);
  markChangedSolutionPK(tag_next_); // mark the solution as changed, as modifying with
                                    // consistent faces will then get the updated boundary
                                    // conditions
  return true;
}

// bool
// Richards::ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u)
// {
//   Teuchos::OSTab tab = vo_->getOSTab();
//   if (vo_->os_OK(Teuchos::VERB_EXTREME))
//     *vo_->os() << "  modifications for consistent face pressures." << std::endl;

//   CalculateConsistentFaces(u->getData().ptr());
//   return true;
// }

// bool
// Richards::ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u)
// {
//   AMANZI_ASSERT(0);
//   return false;
// }


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

// void
// Richards::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u)
// {
//   if (!u->hasComponent("face")) return; // not need

//   // VerboseObject stuff.
//   Teuchos::OSTab tab = vo_->getOSTab();
//   if (vo_->os_OK(Teuchos::VERB_EXTREME))
//     *vo_->os() << "  Modifying predictor for consistent faces" << std::endl;

//   // average cells to faces to give a reasonable initial guess
//   u->scatterMasterToGhosted("cell");
//   const Epetra_MultiVector& u_c = *u->viewComponent("cell", true);
//   Epetra_MultiVector& u_f = *u->viewComponent("face", false);

//   int f_owned = u_f.MyLength();
//   for (int f = 0; f != f_owned; ++f) {
//     AmanziMesh::Entity_ID_List cells;
//     cells = mesh_->getFaceCells(f);
//     int ncells = cells.size();

//     double face_value = 0.0;
//     for (int n = 0; n != ncells; ++n) { face_value += u_c[0][cells[n]]; }
//     u_f[0][f] = face_value / ncells;
//   }
//   markChangedSolution();

//   // Using the old BCs, so should use the old rel perm?
//   // update the rel perm according to the scheme of choice
//   //  UpdatePermeabilityData_(tag_next_);
//   Teuchos::RCP<const CompositeVector> rel_perm =
//     S_->GetPtr<CompositeVector>(uw_coef_key_, tag_next_);

//   S_->GetEvaluator(mass_dens_key_, tag_next_).Update(*S_, name_);
//   Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag_next_);

//   // Update the preconditioner with darcy and gravity fluxes
//   matrix_->Init();
//   matrix_diff_->SetDensity(rho);
//   matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
//   matrix_diff_->UpdateMatrices(Teuchos::null, u);
//   matrix_diff_->ApplyBCs(true, true, true);

//   // derive the consistent faces, involves a solve
//   db_->WriteVector(" p_cf guess:", u.ptr(), true);
//   matrix_diff_->UpdateConsistentFaces(*u);
//   db_->WriteVector(" p_cf soln:", u.ptr(), true);
// }


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

  // // For some reason, wandering PKs break most frequently with an unreasonable
  // // pressure.  This simply tries to catch that before it happens.
  // Teuchos::RCP<const CompositeVector> pres = up->getData();
  // double minT, maxT;

  // const Epetra_MultiVector& pres_c = *pres->viewComponent("cell", false);
  // double minT_c(1.e15), maxT_c(-1.e15);
  // int min_c(-1), max_c(-1);
  // for (int c = 0; c != pres_c.MyLength(); ++c) {
  //   if (pres_c[0][c] < minT_c) {
  //     minT_c = pres_c[0][c];
  //     min_c = c;
  //   }
  //   if (pres_c[0][c] > maxT_c) {
  //     maxT_c = pres_c[0][c];
  //     max_c = c;
  //   }
  // }

  // double minT_f(1.e15), maxT_f(-1.e15);
  // int min_f(-1), max_f(-1);
  // if (pres->hasComponent("face")) {
  //   const Epetra_MultiVector& pres_f = *pres->viewComponent("face", false);
  //   for (int f = 0; f != pres_f.MyLength(); ++f) {
  //     if (pres_f[0][f] < minT_f) {
  //       minT_f = pres_f[0][f];
  //       min_f = f;
  //     }
  //     if (pres_f[0][f] > maxT_f) {
  //       maxT_f = pres_f[0][f];
  //       max_f = f;
  //     }
  //   }
  //   minT = std::min(minT_c, minT_f);
  //   maxT = std::max(maxT_c, maxT_f);

  // } else {
  //   minT = minT_c;
  //   maxT = maxT_c;
  // }

  // double minT_l = minT;
  // double maxT_l = maxT;
  // mesh_->getComm()->MaxAll(&maxT_l, &maxT, 1);
  // mesh_->getComm()->MinAll(&minT_l, &minT, 1);

  // if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //   *vo_->os() << "    Admissible p? (min/max): " << minT << ",  " << maxT << std::endl;
  // }

  // if (minT < -1.e9 || maxT > 1.e8) {
  //   if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
  //     *vo_->os() << " is not admissible, as it is not within bounds of constitutive models:"
  //                << std::endl;

  //     Teuchos::RCP<const Comm_type> comm_p = mesh_->getComm();
  //     Teuchos::RCP<const MpiComm_type> mpi_comm_p =
  //       Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  //     const MPI_Comm& comm = mpi_comm_p->Comm();

  //     ENorm_t global_minT_c, local_minT_c;
  //     ENorm_t global_maxT_c, local_maxT_c;

  //     local_minT_c.value = minT_c;
  //     local_minT_c.gid = pres_c.Map().GID(min_c);
  //     local_maxT_c.value = maxT_c;
  //     local_maxT_c.gid = pres_c.Map().GID(max_c);

  //     MPI_Allreduce(&local_minT_c, &global_minT_c, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
  //     MPI_Allreduce(&local_maxT_c, &global_maxT_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
  //     *vo_->os() << "   cells (min/max): [" << global_minT_c.gid << "] " << global_minT_c.value
  //                << ", [" << global_maxT_c.gid << "] " << global_maxT_c.value << std::endl;

  //     if (pres->hasComponent("face")) {
  //       const Epetra_MultiVector& pres_f = *pres->viewComponent("face", false);
  //       ENorm_t global_minT_f, local_minT_f;
  //       ENorm_t global_maxT_f, local_maxT_f;

  //       local_minT_f.value = minT_f;
  //       local_minT_f.gid = pres_f.Map().GID(min_f);
  //       local_maxT_f.value = maxT_f;
  //       local_maxT_f.gid = pres_f.Map().GID(max_f);

  //       MPI_Allreduce(&local_minT_f, &global_minT_f, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
  //       MPI_Allreduce(&local_maxT_f, &global_maxT_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
  //       *vo_->os() << "   faces (min/max): [" << global_minT_f.gid << "] " << global_minT_f.value
  //                  << ", [" << global_maxT_f.gid << "] " << global_maxT_f.value << std::endl;
  //     }
  //   }
  //   return false;
  // }
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
  if (du->getData()->hasComponent("boundary_face")) {
    du->getData()->getComponent("boundary_face")->putScalar(0.);
  }

  // // debugging -- remove me! --etc
  // for (CompositeVector::name_iterator comp = du->getData()->begin(); comp != du->getData()->end();
  //      ++comp) {
  //   Epetra_MultiVector& du_c = *du->getData()->viewComponent(*comp, false);
  //   double max, l2;
  //   du_c.NormInf(&max);
  //   du_c.Norm2(&l2);
  //   if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //     *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2
  //                << std::endl;
  //   }
  // }

  // // limit by capping corrections when they cross atmospheric pressure
  // // (where pressure derivatives are discontinuous)
  // int my_limited = 0;
  // int n_limited_spurt = 0;
  // if (patm_limit_ > 0.) {
  //   double patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);
  //   for (CompositeVector::name_iterator comp = du->getData()->begin(); comp != du->getData()->end();
  //        ++comp) {
  //     Epetra_MultiVector& du_c = *du->getData()->viewComponent(*comp, false);
  //     const Epetra_MultiVector& u_c = *u->getData()->viewComponent(*comp, false);

  //     for (int c = 0; c != du_c.MyLength(); ++c) {
  //       if ((u_c[0][c] < patm) && (u_c[0][c] - du_c[0][c] > patm + patm_limit_)) {
  //         du_c[0][c] = u_c[0][c] - (patm + patm_limit_);
  //         my_limited++;
  //       } else if ((u_c[0][c] > patm) && (u_c[0][c] - du_c[0][c] < patm - patm_limit_)) {
  //         du_c[0][c] = u_c[0][c] - (patm - patm_limit_);
  //         my_limited++;
  //       }
  //     }
  //   }
  //   mesh_->getComm()->MaxAll(&my_limited, &n_limited_spurt, 1);
  // }

  // if (n_limited_spurt > 0) {
  //   if (vo_->os_OK(Teuchos::VERB_HIGH)) { *vo_->os() << "  limiting the spurt." << std::endl; }
  // }

  // // debugging -- remove me! --etc
  // for (CompositeVector::name_iterator comp = du->getData()->begin(); comp != du->getData()->end();
  //      ++comp) {
  //   Epetra_MultiVector& du_c = *du->getData()->viewComponent(*comp, false);
  //   double max, l2;
  //   du_c.NormInf(&max);
  //   du_c.Norm2(&l2);
  //   if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //     *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2
  //                << std::endl;
  //   }
  // }

  // // Limit based on a max pressure change
  // my_limited = 0;
  // int n_limited_change = 0;
  // if (p_limit_ >= 0.) {
  //   for (CompositeVector::name_iterator comp = du->getData()->begin(); comp != du->getData()->end();
  //        ++comp) {
  //     Epetra_MultiVector& du_c = *du->getData()->viewComponent(*comp, false);

  //     double max;
  //     du_c.NormInf(&max);
  //     if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //       *vo_->os() << "Max pressure correction (" << *comp << ") = " << max << std::endl;
  //     }

  //     for (int c = 0; c != du_c.MyLength(); ++c) {
  //       if (std::abs(du_c[0][c]) > p_limit_) {
  //         du_c[0][c] = ((du_c[0][c] > 0) - (du_c[0][c] < 0)) * p_limit_;
  //         my_limited++;
  //       }
  //     }
  //   }

  //   mesh_->getComm()->MaxAll(&my_limited, &n_limited_change, 1);
  // }

  // // debugging -- remove me! --etc
  // for (CompositeVector::name_iterator comp = du->getData()->begin(); comp != du->getData()->end();
  //      ++comp) {
  //   Epetra_MultiVector& du_c = *du->getData()->viewComponent(*comp, false);
  //   double max, l2;
  //   du_c.NormInf(&max);
  //   du_c.Norm2(&l2);
  //   if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //     *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2
  //                << std::endl;
  //   }
  // }

  // if (n_limited_change > 0) {
  //   if (vo_->os_OK(Teuchos::VERB_HIGH)) { *vo_->os() << "  limited by pressure." << std::endl; }
  // }

  // if (n_limited_spurt > 0) {
  //   return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING;
  // } else if (n_limited_change > 0) {
  //   return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  // }

  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

} // namespace Flow
} // namespace Amanzi
