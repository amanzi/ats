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

#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Accumulation.hh"

#include "upwind_potential_difference.hh"
#include "upwind_cell_centered.hh"
#include "upwind_total_flux.hh"
#include "UpwindFluxFactory.hh"

#include "PK_Helpers.hh"

#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

const std::string OverlandPressureFlow::pk_type_ = "overland flow, pressure basis";


OverlandPressureFlow::OverlandPressureFlow(const Comm_ptr_type& comm,
                                           Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                           const Teuchos::RCP<State>& S)
  : PK_PhysicalBDF_Default(comm, pk_tree, glist, S),
    is_source_term_(false),
    coupled_to_subsurface_via_head_(false),
    coupled_to_subsurface_via_flux_(false),
    perm_update_required_(true),
    precon_used_(true),
    precon_scaled_(false),
    jacobian_(false),
    jacobian_lag_(0),
    iter_(0),
    iter_counter_time_(0.)
{}


void
OverlandPressureFlow::modifyParameterList()
{
  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", 0.01 * 55000.0); // h * nl

  // set some defaults for inherited PKs
  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set<std::string>("conserved quantity key suffix", "water_content");

  PK_PhysicalBDF_Default::modifyParameterList();
}

void
OverlandPressureFlow::parseParameterList()
{
  PK_PhysicalBDF_Default::parseParameterList();

  // get keys
  potential_key_ = Keys::readKey(*plist_, domain_, "potential", "pres_elev");
  flux_key_ = Keys::readKey(*plist_, domain_, "water flux", "water_flux");
  flux_dir_key_ = Keys::readKey(*plist_, domain_, "water flux direction", "water_flux_direction");
  velocity_key_ = Keys::readKey(*plist_, domain_, "velocity", "velocity");
  elev_key_ = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  pd_key_ = Keys::readKey(*plist_, domain_, "ponded depth", "ponded_depth");
  pd_bar_key_ =
    Keys::readKey(*plist_, domain_, "ponded depth, negative", Keys::getVarName(pd_key_) + "_bar");
  wc_bar_key_ = Keys::readKey(
    *plist_, domain_, "conserved quantity, negative", Keys::getVarName(conserved_key_) + "_bar");
  cond_key_ = Keys::readKey(*plist_, domain_, "overland conductivity", "overland_conductivity");
  uw_cond_key_ =
    Keys::readKey(*plist_, domain_, "upwind overland conductivity", "upwind_overland_conductivity");
  molar_dens_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
  mass_dens_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");

  // alter lists for evaluators
  // -- add _bar evaluators
  Teuchos::ParameterList& pd_bar_list = S_->GetEvaluatorList(pd_bar_key_);
  pd_bar_list.setParametersNotAlreadySet(S_->GetEvaluatorList(pd_key_));
  pd_bar_list.set("allow negative ponded depth", true);

  Teuchos::ParameterList& wc_bar_list = S_->GetEvaluatorList(wc_bar_key_);
  wc_bar_list.setParametersNotAlreadySet(S_->GetEvaluatorList(conserved_key_));
  wc_bar_list.set("allow negative water content", true);

  // -- potential evaluator
  auto& potential_list = S_->GetEvaluatorList(potential_key_);
  potential_list.set("evaluator type", "additive");
  potential_list.set<Teuchos::Array<std::string>>("dependencies",
                                                  std::vector<std::string>{ pd_key_, elev_key_ });
  potential_list.set("dependency tags are my tag", true);

  // limiters
  p_limit_ = plist_->get<double>("limit correction to pressure change [Pa]", -1.);
  patm_limit_ =
    plist_->get<double>("limit correction when crossing atmospheric pressure [Pa]", -1.);
  patm_hard_limit_ = plist_->get<bool>("allow no negative ponded depths", false);
  min_vel_ponded_depth_ = plist_->get<double>("min ponded depth for velocity calculation", 1e-2);
  min_tidal_bc_ponded_depth_ = plist_->get<double>("min ponded depth for tidal bc", 0.02);

  // boundary conditions
  auto bc_plist = Teuchos::sublist(plist_, "boundary conditions");
  std::vector<std::string> bcs_found;

  for (auto& bc_sublist : *bc_plist) {
    if (bc_sublist.first == "water level") {
      // -- default Dirichlet conditions
      std::string bc_name = name_ + "_bcs_water_level";
      Teuchos::ParameterList& bc_list = S_->GetEvaluatorList(bc_name);
      bc_list.set("water level", bc_plist->sublist("water level"));
      bc_list.set<std::string>("evaluator type", "independent variable patch")
        .set<std::string>("function list name", "water level")
        .set<std::string>("function inner list name", "boundary water level [m]");
      bcs_found.emplace_back(bc_name);

    } else if (bc_sublist.first == "water flux") {
      // -- default Neumann conditions
      std::string bc_flux_name = name_ + "_bcs_water_flux";
      Teuchos::ParameterList& bc_flux_list = S_->GetEvaluatorList(bc_flux_name);
      bc_flux_list.set("water flux", bc_plist->sublist("water flux"));
      bc_flux_list.set<std::string>("evaluator type", "independent variable patch")
        .set<std::string>("function list name", "water flux")
        .set<std::string>("function inner list name", "outward water flux [mol m^-1 s^-1]");
      bcs_found.emplace_back(bc_flux_name);
    } else if (bc_sublist.first == "ponded depth") {
      // -- Dirichlet with value of this ponded depth + mesh elevation
      auto bc_pd_name = Keys::getKey(domain_, name_ + "_bcs_ponded_depth");
      Teuchos::ParameterList& bc_pd_list = S_->GetEvaluatorList(bc_pd_name);
      bc_pd_list.set("ponded depth", bc_plist->sublist("ponded depth"));
      bc_pd_list.set<std::string>("evaluator type", "flow BC ponded depth");
      bcs_found.emplace_back(bc_pd_name);
    }
    // Need to add:
    //  - seepage face head
    //  - zero gradient
    //  - critical depth
    //  - ???
  }

  // and the boundary condition aggregator
  S_->GetEvaluatorList(name_ + "_bcs")
    .set<std::string>("evaluator type", "boundary condition aggregator")
    .set<std::string>("entity kind", "face")
    .set<Teuchos::Array<std::string>>("dependencies", bcs_found);
}


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void
OverlandPressureFlow::setup()
{
  PK_PhysicalBDF_Default::setup();

  // -- water content, and evaluator, and derivative for PC
  PKHelpers::requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  //    and at the current time, where it is a copy evaluator
  PKHelpers::requireAtCurrent(conserved_key_, tag_current_, *S_, name_);

  // this pk uses density to invert for velocity from flux
  PKHelpers::requireAtNext(molar_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  SetupOverlandFlow_();
  SetupPhysicalEvaluators_();
}


void
OverlandPressureFlow::SetupOverlandFlow_()
{
  //
  // Diffusion Operators
  // ------------------------------------------------------------------
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  auto& bc_fac = S_->Require<Operators::BCs, Operators::BCs_Factory>(name_ + "_bcs", tag_next_);
  bc_fac.set_mesh(mesh_);
  bc_fac.set_entity_kind(AmanziMesh::Entity_kind::FACE);
  bc_fac.set_dof_type(WhetStone::DOF_Type::SCALAR);
  S_->RequireEvaluator(name_ + "_bcs", tag_next_);

  if (plist_->sublist("boundary conditions").isSublist("water level"))
    S_->Require<MultiPatch<double>, MultiPatchSpace>(name_ + "_bcs_water_level", tag_next_)
      .set_flag(Operators::OPERATOR_BC_DIRICHLET);
  if (plist_->sublist("boundary conditions").isSublist("water flux"))
    S_->Require<MultiPatch<double>, MultiPatchSpace>(name_ + "_bcs_water_flux", tag_next_)
      .set_flag(Operators::OPERATOR_BC_NEUMANN);

  // bcs require this on boundary_faces for pd
  // S_->Require<CompositeVector, CompositeVectorSpace>(elev_key_, tag_next_)
  //   .SetMesh(mesh_)
  //   ->SetGhosted(true)
  //   ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
  //   ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  // S_->RequireEvaluator(elev_key_, tag_next_);

  // -- nonlinear coefficients and upwinding
  Teuchos::ParameterList& upwind_plist = plist_->sublist("upwinding");
  upwinding_ = Operators::UpwindFactory::Create(upwind_plist, *S_, name_, tag_next_, flux_dir_key_);

  std::string coef_location = upwinding_->CoefficientLocation();
  PKHelpers::requireNonlinearDiffusionCoefficient(uw_cond_key_, tag_next_, coef_location, *S_);

  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  if (!mfd_plist.isParameter("scaled constraint equation"))
    mfd_plist.set("scaled constraint equation", true);

  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.Create(mfd_plist, mesh_);
  matrix_diff_->SetTensorCoefficient(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // -- create the operator, data for flux directions
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = opfactory.Create(face_diff_list, mesh_, bc_);
  face_matrix_diff_->SetTensorCoefficient(Teuchos::null);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  face_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  S_->Require<CompositeVector, CompositeVectorSpace>(flux_dir_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // -- create the operators for the preconditioner
  //    diffusion
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  mfd_pc_plist.set("scaled constraint equation", mfd_plist.get<bool>("scaled constraint equation"));
  mfd_pc_plist.set("constraint equation scaling cutoff",
                   mfd_plist.get<double>("constraint equation scaling cutoff", 1.0));
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

  //    inverse for preconditioner
  precon_used_ = plist_->isSublist("inverse");
  if (precon_used_) {
    auto& inv_list = mfd_pc_plist.sublist("inverse");
    inv_list.setParameters(plist_->sublist("inverse"));
  }
  precon_scaled_ = plist_->get<bool>("scale preconditioner to pressure", !precon_used_);

  //    create the preconditioner -- can this be done via clone?
  preconditioner_diff_ = opfactory.Create(mfd_pc_plist, mesh_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  // -- if using approximate Jacobian for the preconditioner, we also need
  // -- derivative information.
  jacobian_ = (mfd_pc_plist.get<std::string>("Newton correction", "none") != "none");
  if (jacobian_) {
    jacobian_lag_ = mfd_pc_plist.get<int>("Newton correction lag", 0);

    // require the derivative drel_perm/dp
    S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
        cond_key_, tag_next_, pd_key_, tag_next_)
      .SetGhosted();
    if (mfd_pc_plist.get<std::string>("discretization primary") != "fv: default") {
      // MFD -- upwind required, require data
      duw_cond_key_ = Keys::getDerivKey(uw_cond_key_, pd_key_);

      // note this is always done on faces using total flux upwinding
      PKHelpers::requireNonlinearDiffusionCoefficient(
        duw_cond_key_, tag_next_, "upwind: face", *S_);

      double flux_eps = plist_->get<double>("upwind flux epsilon", 1.e-12);
      upwinding_dkdp_ =
        Teuchos::rcp(new Operators::UpwindTotalFlux(name_, tag_next_, flux_dir_key_, flux_eps));
    } else {
      // FV -- no upwinding of derivative
      duw_cond_key_ = std::string();
    }
  }

  //
  // Accumluation Operator
  // ------------------------------------------------------------------
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ =
    Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  //
  // Source term
  // ------------------------------------------------------------------
  is_source_term_ = plist_->get<bool>("source term", false);
  if (is_source_term_) {
    if (source_key_.empty()) {
      source_key_ = Keys::readKey(*plist_, domain_, "source", "water_source");
    }

    bool source_in_meters = plist_->get<bool>("water source in meters", true);
    if (source_in_meters) {
      // create multiplicative evaluator for the product of density and the source
      Key source_key_meters = source_key_;
      source_key_ = source_key_meters + "_mols_per_s";

      S_->GetEvaluatorList(source_key_)
        .set<std::string>("evaluator type", "multiplicative")
        .set<Teuchos::Array<std::string>>(
          "dependencies", std::vector<std::string>{ source_key_meters, molar_dens_key_ });
    }

    PKHelpers::requireAtNext(source_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  //
  // Coupling to Subsurface
  // ------------------------------------------------------------------
  coupled_to_subsurface_via_flux_ = plist_->get<bool>("coupled to subsurface via flux", false);
  coupled_to_subsurface_via_head_ = plist_->get<bool>("coupled to subsurface via head", false);
  AMANZI_ASSERT(!(coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_head_));

  if (coupled_to_subsurface_via_head_) {
    // -- source term from subsurface, filled in by evaluator,
    //    which picks the fluxes from "water_flux" field.
    ss_flux_key_ = Keys::readKey(*plist_, domain_, "surface-subsurface flux", "subsurface_flux");
    PKHelpers::requireAtNext(ss_flux_key_, tag_next_, *S_)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  //
  // Require fields and evaluators
  // ------------------------------------------------------------------
  //  NOTE: no need to require evaluator for p here, this was done in pk_physical
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_)
    .Update(*matrix_->getRangeMap())
    ->SetGhosted()
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

  // potential may not actually need cells, but for debugging and sanity's sake, we require them
  PKHelpers::requireAtNext(potential_key_, tag_next_, *S_)
    .Update(*matrix_->getRangeMap())
    ->SetGhosted()
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

  // flux
  PKHelpers::requireAtNext(flux_key_, tag_next_, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  // velocity for diagnostics
  PKHelpers::requireAtNext(velocity_key_, Tags::NEXT, *S_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 3);
};


// -----------------------------------------------------------------------
// Create the physical evaluators for things specific to overland_pressure
// -----------------------------------------------------------------------
void
OverlandPressureFlow::SetupPhysicalEvaluators_()
{
  // -- water content bar (can be negative)
  PKHelpers::requireAtNext(wc_bar_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    wc_bar_key_, tag_next_, key_, tag_next_);

  // -- ponded depth
  PKHelpers::requireAtNext(pd_key_, tag_next_, *S_).Update(*matrix_->getRangeMap())->SetGhosted();
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(pd_key_, tag_next_, key_, tag_next_);
  //    ...with a copy at the old time
  PKHelpers::requireAtCurrent(pd_key_, tag_current_, *S_, pd_key_);

  // -- ponded depth bar (can be negative)
  PKHelpers::requireAtNext(pd_bar_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    pd_bar_key_, tag_next_, key_, tag_next_);

  // -- conductivity evaluator
  PKHelpers::requireAtNext(cond_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void
OverlandPressureFlow::initialize()
{
  // Initialize BDF stuff and physical domain stuff.
  PK_PhysicalBDF_Default::initialize();

  // Set extra fields as initialized -- these don't currently have evaluators.
  S_->GetW<CompositeVector>(uw_cond_key_, tag_next_, uw_cond_key_).putScalar(0.0);
  S_->GetRecordW(uw_cond_key_, tag_next_, uw_cond_key_).set_initialized();

  if (!duw_cond_key_.empty()) {
    S_->GetW<CompositeVector>(duw_cond_key_, tag_next_, name_).putScalar(1.0);
    S_->GetRecordW(duw_cond_key_, tag_next_, name_).set_initialized();
  }

  S_->GetW<CompositeVector>(flux_key_, tag_next_, getName()).putScalar(0.0);
  S_->GetRecordW(flux_key_, tag_next_, getName()).set_initialized();
  PKHelpers::changedEvaluatorPrimary(flux_key_, tag_next_, *S_);

  S_->GetW<CompositeVector>(flux_dir_key_, tag_next_, name_).putScalar(0.);
  S_->GetRecordW(flux_dir_key_, tag_next_, name_).set_initialized();

  S_->GetW<CompositeVector>(velocity_key_, Tags::NEXT, name_).putScalar(0.);
  PKHelpers::changedEvaluatorPrimary(velocity_key_, Tags::NEXT, *S_);
  S_->GetRecordW(velocity_key_, Tags::NEXT, name_).set_initialized();

  // operators
  const auto& bc = S_->GetPtr<Operators::BCs>(name_ + "_bcs", tag_next_);
  matrix_diff_->SetBCs(bc, bc);
  preconditioner_diff_->SetBCs(bc, bc);
  face_matrix_diff_->SetBCs(bc, bc);
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void
OverlandPressureFlow::commitStep(double t_old, double t_new, const Tag& tag_next)
{
  // saves primary variable, conserved quantity
  PK_PhysicalBDF_Default::commitStep(t_old, t_new, tag_next);

  // also save ponded depth
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(pd_key_, tag_current, tag_next, *S_);
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool
OverlandPressureFlow::UpdatePermeabilityData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability?";

  bool update_perm = S_->GetEvaluator(pd_key_, tag).Update(*S_, name_);
  update_perm |= S_->GetEvaluator(potential_key_, tag).Update(*S_, name_);
  update_perm |= S_->GetEvaluator(cond_key_, tag).Update(*S_, name_);
  update_perm |= perm_update_required_;

  if (update_perm) {
    // Update the perm only if needed.
    perm_update_required_ = false;

    // update the direction of the flux -- note this is NOT the flux, but grad(h+z).
    Teuchos::RCP<CompositeVector> flux_dir =
      S_->GetPtrW<CompositeVector>(flux_dir_key_, tag, name_);
    Teuchos::RCP<const CompositeVector> pres_elev =
      S_->GetPtr<CompositeVector>(potential_key_, tag);

    // -- need to apply BCs to get boundary flux directions correct
    FixBCsForOperator_(tag, face_matrix_diff_.ptr()); // deals with zero gradient condition

    // -- now we can calculate the flux direction
    // face_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null); // can this be removed?
    face_matrix_diff_->UpdateFlux(pres_elev.ptr(), flux_dir.ptr());
    flux_dir->scatterMasterToGhosted("face"); // for vis only
    //face_matrix_diff_->ApplyBCs(true, true, true); // because we don't do this....

    // upwind
    // -- get upwind conductivity data
    Teuchos::RCP<CompositeVector> uw_cond =
      S_->GetPtrW<CompositeVector>(uw_cond_key_, tag, uw_cond_key_);

    // -- get conductivity data
    Teuchos::RCP<const CompositeVector> cond = S_->GetPtr<CompositeVector>(cond_key_, tag);

    // -- Move rel perm on boundary_faces into uw_rel_perm on faces
    // Move rel perm on boundary_faces into uw_rel_perm on faces
    const auto& vandelay = mesh_->getBoundaryFaceImporter();
    {
      const auto& cond_bf = *cond->getComponent("boundary_face", false);
      auto& uw_cond_f = *uw_cond->getComponent("face", false);
      uw_cond_f.doExport(cond_bf, vandelay, Tpetra::CombineMode::INSERT);
    }

    // -- upwind
    upwinding_->Update(*cond, *uw_cond, *S_);
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Derivatives of the overland conductivity, upwinded.
// -----------------------------------------------------------------------------
bool
OverlandPressureFlow::UpdatePermeabilityDerivativeData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability derivatives?";

  bool update_perm = S_->GetEvaluator(cond_key_, tag).UpdateDerivative(*S_, name_, pd_key_, tag);
  Teuchos::RCP<const CompositeVector> dcond =
    S_->GetDerivativePtr<CompositeVector>(cond_key_, tag, pd_key_, tag);

  if (update_perm) {
    if (!duw_cond_key_.empty()) {
      // get upwind conductivity data
      Teuchos::RCP<CompositeVector> duw_cond =
        S_->GetPtrW<CompositeVector>(duw_cond_key_, tag, name_);
      duw_cond->putScalar(0.);

      // Then upwind.  This overwrites the boundary if upwinding says so.
      upwinding_dkdp_->Update(*dcond, *duw_cond, *S_);
    }
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the tag
// -----------------------------------------------------------------------------
void
OverlandPressureFlow::UpdateBoundaryConditions_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating BCs." << std::endl;

  S_->GetEvaluator(name_ + "_bcs", tag_next_).Update(*S_, name_);
}


// -----------------------------------------------------------------------------
// Given a vector, apply the Dirichlet data to that vector.
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void
OverlandPressureFlow::ApplyDirichletBCs_(const Operators::BCs& bcs,
                                         CompositeVector& u,
                                         const CompositeVector& elev)
{
  if (u.hasComponent("face")) {
    auto u_f = u.viewComponent("face", false);
    auto elev_f = u.viewComponent("face", false);
    auto bc_value = bcs.bc_value();
    auto bc_model = bcs.bc_model();
    Kokkos::parallel_for(
      "applyDirichletBCs", u_f.extent(0), KOKKOS_LAMBDA(const int& f) {
        if (bc_model(f) == Operators::OPERATOR_BC_DIRICHLET) {
          u_f(f, 0) = bc_value(f) - elev_f(f, 0);
        }
      });
  }

  if (u.hasComponent("boundary_face")) {
    auto u_f = u.viewComponent("boundary_face", false);
    auto elev_f = u.viewComponent("face", false);
    auto bc_value = bcs.bc_value();
    auto bc_model = bcs.bc_model();

    const AmanziMesh::MeshCache& mesh = u.getMesh()->getCache();
    Kokkos::parallel_for(
      "applyDirichletBCs", u_f.extent(0), KOKKOS_LAMBDA(const int& bf) {
        auto f = mesh.getBoundaryFaceFace(bf);
        if (bc_model(f) == Operators::OPERATOR_BC_DIRICHLET) {
          u_f(bf, 0) = bc_value(f) - elev_f(f, 0);
        }
      });
  }
}


void
OverlandPressureFlow::FixBCsForOperator_(const Tag& tag,
                                         const Teuchos::Ptr<Operators::PDE_Diffusion>& diff_op)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the Operator." << std::endl;

  // auto& markers = bc_markers();
  // auto& values = bc_values();

  // // Now we can safely calculate q = -k grad z for zero-gradient problems
  // if (bc_zero_gradient_->size() > 0) {
  //   Teuchos::RCP<const CompositeVector> elev = S_->GetPtr<CompositeVector>(elev_key_, tag);
  //   Teuchos::RCP<const CompositeVector> ponded_depth = S_->GetPtr<CompositeVector>(pd_key_, tag);
  //   Teuchos::RCP<const CompositeVector> upw_conductivity =
  //     S_->GetPtr<CompositeVector>(uw_cond_key_, tag);
  //   Teuchos::RCP<const CompositeVector> conductivity = S_->GetPtr<CompositeVector>(cond_key_, tag);
  //   Teuchos::RCP<const CompositeVector> flux = S_->GetPtr<CompositeVector>(flux_key_, tag);

  //   const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face", false);
  //   const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell", false);
  //   const Epetra_MultiVector& ponded_c = *ponded_depth->ViewComponent("cell", false);
  //   const Epetra_MultiVector& upw_cond_f = *upw_conductivity->ViewComponent("face", false);
  //   const Epetra_MultiVector& cond_c = *conductivity->ViewComponent("cell", false);
  //   const Epetra_MultiVector& flux_f = *flux->ViewComponent("face", false);

  //   std::vector<WhetStone::DenseMatrix>& Aff = diff_op->local_op()->matrices;

  //   int ncells_owned = mesh_->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  //   int nfaces_owned = mesh_->num_entities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  //   for (const auto& bc : *bc_zero_gradient_) {
  //     int f = bc.first;

  //     AmanziMesh::Entity_ID_List cells;
  //     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  //     AMANZI_ASSERT(cells.size() == 1);
  //     AmanziMesh::Entity_ID c = cells[0];

  //     if (f < nfaces_owned) {
  //       double dp = elevation_f[0][f] - elevation_c[0][c];
  //       double bc_val = -dp * Aff[f](0, 0);

  //       markers[f] = Operators::OPERATOR_BC_NEUMANN;
  //       values[f] = bc_val / mesh_->face_area(f);
  //     }
  //   }
  // }

  // if (bc_tidal_->size() > 0) {
  //   const double& p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  //   const Epetra_MultiVector& rho_l =
  //     *S_->Get<CompositeVector>(mass_dens_key_, tag).ViewComponent("cell");
  //   Teuchos::RCP<const CompositeVector> elev = S_->GetPtr<CompositeVector>(elev_key_, tag);
  //   const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face", false);
  //   const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell", false);

  //   std::vector<WhetStone::DenseMatrix>& Aff = diff_op->local_op()->matrices;
  //   double gz = -(S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT))[2];
  //   int nfaces_owned = mesh_->num_entities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);

  //   for (const auto& bc : *bc_tidal_) {
  //     int f = bc.first;

  //     AmanziMesh::Entity_ID_List cells, faces;
  //     std::vector<int> fdirs;
  //     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  //     AMANZI_ASSERT(cells.size() == 1);
  //     AmanziMesh::Entity_ID c = cells[0];

  //     if (f < nfaces_owned) {
  //       mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  //       int j = 0;
  //       for (j = 0; j < faces.size(); j++) {
  //         if (faces[j] == f) break;
  //       }
  //       if (j >= faces.size()) {
  //         Errors::Message msg(
  //           "Overland_pressure PK: boundary face is not found in the boundary cell.");
  //         Exceptions::amanzi_throw(msg);
  //       }
  //       double h0 = bc.second;

  //       if ((h0 - elevation_f[0][f] < min_tidal_bc_ponded_depth_)) {
  //         double dp = elevation_f[0][f] - elevation_c[0][c];
  //         double bc_val = -dp * Aff[f](0, 0);

  //         markers[f] = Operators::OPERATOR_BC_NEUMANN;
  //         values[f] = bc_val / mesh_->face_area(f);
  //       } else {
  //         markers[f] = Operators::OPERATOR_BC_DIRICHLET;
  //         values[f] = h0;
  //       }

  //       if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
  //         *vo_->os() << "Tidal BC2: f=" << f << " type " << markers[f] << " val " << values[f]
  //                    << " Aff " << Aff[f](0, 0) << "\n";
  //       }
  //     }
  //   }
  // }
};


void
OverlandPressureFlow::FixBCsForPrecon_(const Tag& tag)
{}

bool
OverlandPressureFlow::ModifyPredictor(double h,
                                      Teuchos::RCP<const TreeVector> u0,
                                      Teuchos::RCP<TreeVector> u)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Modifying predictor:" << std::endl;

  // update boundary conditions
  UpdateBoundaryConditions_(tag_next_);
  // db_->WriteBoundaryConditions(bc_->model(), bc_->value());

  // push Dirichlet data into predictor
  const auto& bcs = S_->Get<Operators::BCs>(name_ + "_bcs", tag_next_);
  PKHelpers::applyDirichletBCs(bcs, *u->getData());

  return false;
};


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
OverlandPressureFlow::ModifyCorrection(double h,
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

  // limit by capping corrections when they cross atmospheric pressure
  // (where pressure derivatives are discontinuous)
  int my_limited = 0;
  int n_limited_spurt = 0;
  if (patm_limit_ > 0.) {
    double patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    auto du_c = du->getData()->viewComponent("cell", false);
    const auto u_c = u->getData()->viewComponent("cell", false);
    double patm_limit(patm_limit_);

    Kokkos::parallel_reduce(
      "OverlandFlowPressure::ModifyPredictor 'limit correction when crossing atmospheric pressure "
      "[Pa]'",
      du_c.extent(0),
      KOKKOS_LAMBDA(const int& c, int& count) {
        if ((u_c(c, 0) < patm) && (u_c(c, 0) - du_c(c, 0) > patm + patm_limit)) {
          du_c(c, 0) = u_c(c, 0) - (patm + patm_limit);
          count++;
        } else if ((u_c(c, 0) > patm) && (u_c(c, 0) - du_c(c, 0) < patm - patm_limit)) {
          du_c(c, 0) = u_c(c, 0) - (patm - patm_limit);
          count++;
        }
      },
      my_limited);
    Teuchos::reduceAll(*mesh_->getComm(), Teuchos::REDUCE_SUM, 1, &my_limited, &n_limited_spurt);
  }

  if (patm_hard_limit_) {
    double patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    auto du_c = du->getData()->viewComponent("cell", false);
    const auto u_c = u->getData()->viewComponent("cell", false);

    Kokkos::parallel_reduce(
      "OverlandFlowPressure::ModifyPredictor 'allow no negative ponded depths'",
      du_c.extent(0),
      KOKKOS_LAMBDA(const int& c, int& count) {
        if (u_c(c, 0) - du_c(c, 0) < patm) {
          du_c(c, 0) = u_c(c, 0) - patm;
          count++;
        }
      },
      my_limited);
    Teuchos::reduceAll(*mesh_->getComm(), Teuchos::REDUCE_SUM, 1, &my_limited, &n_limited_spurt);
  }

  if (n_limited_spurt > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) { *vo_->os() << "  limiting the spurt." << std::endl; }
  }

  // Limit based on a max pressure change
  my_limited = 0;
  int n_limited_change = 0;

  if (p_limit_ > 0.) {
    for (const auto& comp : *du->getData()) {
      auto du_c = du->getData()->viewComponent("cell", false);
      const auto u_c = u->getData()->viewComponent("cell", false);
      double p_limit(p_limit_);

      Kokkos::parallel_reduce(
        "OverlandFlowPressure::ModifyPredictor 'allow no negative ponded depths'",
        du_c.extent(0),
        KOKKOS_LAMBDA(const int& c, int& count) {
          if (fabs(du_c(c, 0)) > p_limit) {
            du_c(c, 0) = ((du_c(c, 0) > 0) - (du_c(c, 0) < 0)) * p_limit;
            count++;
          }
        },
        my_limited);
    }

    Teuchos::reduceAll(*mesh_->getComm(), Teuchos::REDUCE_SUM, 1, &my_limited, &n_limited_change);
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
