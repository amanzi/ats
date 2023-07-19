/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "flow_bc_factory.hh"
#include "Mesh.hh"
#include "Point.hh"
#include "Op.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"
#include "independent_variable_field_evaluator.hh"

#include "upwind_potential_difference.hh"
#include "upwind_total_flux.hh"
#include "upwind_cell_centered.hh"
#include "pres_elev_evaluator.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"

#include "UpwindFluxFactory.hh"

#include "PDE_DiffusionFactory.hh"
#include "overland.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_RES_FLAG 0

OverlandFlow::OverlandFlow(Teuchos::ParameterList& FElist,
                           const Teuchos::RCP<Teuchos::ParameterList>& plist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution)
  : PK(FElist, plist, S, solution),
    PK_PhysicalBDF_Default(FElist, plist, S, solution),
    standalone_mode_(false),
    is_source_term_(false),
    niter_(0)
{
  // used for error norm
  if (!plist_->isParameter("conserved quantity key suffix"))
    plist_->set("conserved quantity key suffix", "ponded_depth_times_cell_volume");

  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .01); // h
}


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void
OverlandFlow::Setup(const Teuchos::Ptr<State>& S)
{
  // set up the meshes
  if (!S->HasMesh("surface")) {
    Teuchos::RCP<const AmanziMesh::Mesh> domain = S->GetMesh();
    //    AMANZI_ASSERT(domain->getSpaceDimension() == 2);
    standalone_mode_ = true;
    S->AliasMesh("domain", "surface");
  } else {
    standalone_mode_ = false;
  }

  PK_PhysicalBDF_Default::Setup(S);
  SetupOverlandFlow_(S);
  SetupPhysicalEvaluators_(S);
}


void
OverlandFlow::SetupOverlandFlow_(const Teuchos::Ptr<State>& S)
{
  // Default evaluators
  S->RequireEvaluator("surface-cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure", Tags::DEFAULT);

  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_head_ = bc_factory.CreateHead();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();
  bc_seepage_head_ = bc_factory.CreateWithFunction("seepage face head", "boundary head");
  bc_critical_depth_ = bc_factory.CreateCriticalDepth();

  // -- nonlinear coefficients/upwinding
  Teuchos::ParameterList& cond_plist = plist_->sublist("overland conductivity evaluator");
  Operators::UpwindFluxFactory upwfactory;

  upwinding_ = upwfactory.Create(cond_plist,
                                 S,
                                 name_,
                                 Keys::getKey(domain_, "overland_conductivity"),
                                 Keys::getKey(domain_, "upwind_overland_conductivity"),
                                 "surface-water_flux_direction");

  // -- require the data on appropriate locations
  std::string coef_location = upwinding_->CoefficientLocation();
  if (coef_location == "upwind: face") {
    S->Require<CompositeVector, CompositeVectorSpace>(
       Keys::getKey(domain_, Tags::NEXT, "upwind_overland_conductivity"), name_)
      ->SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S->Require<CompositeVector, CompositeVectorSpace>(
       Keys::getKey(domain_, Tags::NEXT, "upwind_overland_conductivity"), name_)
      ->SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  } else {
    Errors::Message message("Unknown upwind coefficient location in overland flow.");
    Exceptions::amanzi_throw(message);
  }
  S->GetField(Keys::getKey(domain_, "upwind_overland_conductivity"), name_)->set_io_vis(false);

  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  if (!mfd_plist.isParameter("scaled constraint equation"))
    mfd_plist.set("scaled constraint equation", true);

  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.Create(mfd_plist, mesh_, bc_);
  matrix_diff_->SetTensorCoefficient(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // -- create the operator, data for flux directions
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = opfactory.Create(face_diff_list, mesh_, bc_);
  face_matrix_diff_->SetTensorCoefficient(Teuchos::null);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  face_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  S->Require<CompositeVector, CompositeVectorSpace>(
     "surface-water_flux_direction", Tags::NEXT, name_)
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

  mfd_pc_plist.set("inverse", plist_->sublist("inverse"));
  // old style... deprecate me!
  mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("preconditioner"));
  mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("linear solver"));

  preconditioner_diff_ = opfactory.Create(mfd_pc_plist, mesh_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  // If using approximate Jacobian for the preconditioner, we also need derivative information.
  jacobian_ = mfd_pc_plist.get<std::string>("Newton correction", "none") != "none";
  if (jacobian_) {
    if (preconditioner_->RangeMap().HasComponent("face")) {
      // MFD -- upwind required
      Key duwkey = Keys::getDerivKey(Keys::getKey(domain_, "upwind_overland_conductivity"), key_);
      S->Require<CompositeVector, CompositeVectorSpace>(duwkey, Tags::NEXT, name_)
        ->SetMesh(mesh_)
        ->SetGhosted()
        ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

      Key dkey = Keys::getDerivKey(Keys::getKey(domain_, "overland_conductivity"), key_);
      upwinding_dkdp_ = Teuchos::rcp(
        new Operators::UpwindTotalFlux(name_, dkey, duwkey, "surface-water_flux_direction", 1.e-8));
    }
  }

  //    accumulation
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ =
    Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  // symbolic structure is set

  // primary variable
  S->Require<CompositeVector, CompositeVectorSpace>(
     Keys::getKey(domain_, Tags::NEXT, "ponded_depth"), name_)
    ->Update(matrix_->RangeMap())
    ->SetGhosted();

  // fluxes
  S->Require<CompositeVector, CompositeVectorSpace>("surface-water_flux", Tags::NEXT, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void
OverlandFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S)
{
  // -- evaluator for surface geometry.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2, 1);
  locations2[0] = AmanziMesh::Entity_kind::CELL;
  locations2[1] = AmanziMesh::Entity_kind::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->Require<CompositeVector, CompositeVectorSpace>(Keys::getKey(domain_, Tags::NEXT, "elevation"))
    ->SetMesh(S->GetMesh("surface"))
    ->SetGhosted()
    ->AddComponents(names2, locations2, num_dofs2);

  S->Require<CompositeVector, CompositeVectorSpace>(
     Keys::getKey(domain_, Tags::NEXT, "slope_magnitude"))
    ->SetMesh(S->GetMesh("surface"))
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  Teuchos::RCP<Flow::ElevationEvaluator> elev_evaluator;
  if (standalone_mode_) {
    AMANZI_ASSERT(plist_->isSublist("elevation evaluator"));
    Teuchos::ParameterList elev_plist = plist_->sublist("elevation evaluator");
    elev_plist.set("evaluator name", Keys::getKey(domain_, "elevation"));
    elev_evaluator = Teuchos::rcp(new Flow::StandaloneElevationEvaluator(elev_plist));
  } else {
    Teuchos::ParameterList elev_plist = plist_->sublist("elevation evaluator");
    elev_plist.set("evaluator name", Keys::getKey(domain_, "elevation"));
    elev_evaluator = Teuchos::rcp(new Flow::MeshedElevationEvaluator(elev_plist));
  }
  S->SetEvaluator(Keys::getKey(domain_, "elevation"), elev_evaluator);
  S->SetEvaluator(Keys::getKey(domain_, "slope_magnitude"), elev_evaluator);

  // -- evaluator for potential field, h + z
  S->Require<CompositeVector, CompositeVectorSpace>(Keys::getKey(domain_, Tags::NEXT, "pres_elev"))
    ->Update(matrix_->RangeMap())
    ->SetGhosted();
  Teuchos::ParameterList pres_elev_plist = plist_->sublist("potential evaluator");
  pres_elev_plist.set("evaluator name", Keys::getKey(domain_, "pres_elev"));
  Teuchos::RCP<Flow::PresElevEvaluator> pres_elev_eval =
    Teuchos::rcp(new Flow::PresElevEvaluator(pres_elev_plist));
  S->SetEvaluator(Keys::getKey(domain_, "pres_elev"), pres_elev_eval);

  // -- evaluator for source term
  is_source_term_ = plist_->get<bool>("source term");
  if (is_source_term_) {
    // source term itself [m/s]
    source_key_ = Keys::readKey(*plist_, domain_, "source", "water_source");
    S->Require<CompositeVector, CompositeVectorSpace>(source_key_, Tags::NEXT)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S->RequireEvaluator(source_key_);
  }

  // -- conductivity evaluator
  S->Require<CompositeVector, CompositeVectorSpace>(
     Keys::getKey(domain_, Tags::NEXT, "overland_conductivity"))
    ->SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  AMANZI_ASSERT(plist_->isSublist("overland conductivity evaluator"));
  Teuchos::ParameterList cond_plist = plist_->sublist("overland conductivity evaluator");
  cond_plist.set("evaluator name", Keys::getKey(domain_, "overland_conductivity"));

  Teuchos::RCP<Flow::OverlandConductivityEvaluator> cond_evaluator =
    Teuchos::rcp(new Flow::OverlandConductivityEvaluator(cond_plist));
  S->SetEvaluator(Keys::getKey(domain_, "overland_conductivity"), cond_evaluator);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void
OverlandFlow::Initialize(const Teuchos::Ptr<State>& S)
{
  // Initialize BDF stuff and physical domain stuff.
  PK_PhysicalBDF_Default::Initialize(S);

  // Initialize BC values
  bc_head_->Compute(S->get_time());
  bc_zero_gradient_->Compute(S->get_time());
  bc_flux_->Compute(S->get_time());
  bc_seepage_head_->Compute(S->get_time());
  bc_critical_depth_->Compute(S->get_time());

  // Set extra fields as initialized -- these don't currently have evaluators.
  Key uwkey = Keys::getKey(domain_, "upwind_overland_conductivity");
  S->GetW<CompositeVector>(uwkey, name_).PutScalar(1.0);
  S->GetField(uwkey, name_)->set_initialized();
  if (jacobian_ && preconditioner_->RangeMap().HasComponent("face")) {
    Key dkey = Keys::getDerivKey(uwkey, key_);
    S->GetW<CompositeVector>(dkey, name_).PutScalar(1.0);
    S->GetField(dkey, name_)->set_initialized();
  }
  S->GetField("surface-water_flux", name_)->set_initialized();
  S->GetW<CompositeVector>("surface-water_flux_direction", name_).PutScalar(0.);
  S->GetField("surface-water_flux_direction", name_)->set_initialized();
  //  S->GetField("surface-velocity", name_)->set_initialized();
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void
OverlandFlow::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  niter_ = 0;
  double dt = t_new - t_old;
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Commiting state." << std::endl;

  //PKPhysicalBDFBase::commit_state(, S);
  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, S);

  // update boundary conditions
  bc_head_->Compute(S->get_time());
  bc_flux_->Compute(S->get_time());
  bc_seepage_head_->Compute(S->get_time());
  bc_critical_depth_->Compute(S->get_time());
  UpdateBoundaryConditions_(S.ptr());

  // Update flux if rel perm or h + Z has changed.
  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetEvaluator(Keys::getKey(domain_, "pres_elev"))->HasFieldChanged(S.ptr(), name_);

  // update the stiffness matrix with the new rel perm
  Teuchos::RCP<const CompositeVector> conductivity =
    S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "upwind_overland_conductivity"));
  matrix_->Init();
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // Patch up BCs for zero-gradient
  FixBCsForOperator_(S.ptr());

  // derive the fluxes
  Teuchos::RCP<const CompositeVector> potential =
    S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "pres_elev"));
  Teuchos::RCP<CompositeVector> flux = S->GetPtrW<CompositeVector>("surface-water_flux", name_);
  matrix_diff_->UpdateFlux(potential.ptr(), flux.ptr());
};


// -----------------------------------------------------------------------------
// Update diagnostics -- used prior to vis.
// -----------------------------------------------------------------------------
void
OverlandFlow::CalculateDiagnostics(const Teuchos::RCP<State>& S){};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool
OverlandFlow::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability?";

  bool update_perm =
    S->GetEvaluator(Keys::getKey(domain_, "overland_conductivity"))->HasFieldChanged(S, name_);
  update_perm |= S->GetEvaluator(Keys::getKey(domain_, "ponded_depth"))->HasFieldChanged(S, name_);
  update_perm |= S->GetEvaluator(Keys::getKey(domain_, "pres_elev"))->HasFieldChanged(S, name_);

  if (update_perm) {
    // get upwind conductivity data
    Teuchos::RCP<CompositeVector> uw_cond =
      S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "upwind_overland_conductivity"), name_);

    // update the direction of the flux -- note this is NOT the flux
    Teuchos::RCP<CompositeVector> flux_dir =
      S->GetPtrW<CompositeVector>("surface-water_flux_direction", name_);
    Teuchos::RCP<const CompositeVector> pres_elev =
      S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "pres_elev"));
    face_matrix_diff_->UpdateFlux(pres_elev.ptr(), flux_dir.ptr());

    // get conductivity data
    Teuchos::RCP<const CompositeVector> cond =
      S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "overland_conductivity"));
    const Epetra_MultiVector& cond_c = *cond->ViewComponent("cell", false);

    // place internal cell's value on faces -- this should be fixed to be the boundary data
    { // place boundary_faces on faces
      Epetra_MultiVector& uw_cond_f = *uw_cond->ViewComponent("face", false);

      AmanziMesh::Entity_ID_List cells;
      int nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
      for (int f = 0; f != nfaces_owned; ++f) {
        cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        if (cells.size() == 1) {
          int c = cells[0];
          uw_cond_f[0][f] = cond_c[0][c];
        }
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);
    if (uw_cond->HasComponent("face")) uw_cond->ScatterMasterToGhosted("face");
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Derivatives of the overland conductivity, upwinded.
// -----------------------------------------------------------------------------
bool
OverlandFlow::UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability derivatives?";

  bool update_perm = S->GetEvaluator(Keys::getKey(domain_, "overland_conductivity"))
                       ->HasFieldDerivativeChanged(S, name_, key_);
  Key dcond_key = Keys::getDerivKey(Keys::getKey(domain_, "overland_conductivity"), key_);
  Teuchos::RCP<const CompositeVector> dcond = S->GetPtr<CompositeVector>(dcond_key);

  if (update_perm) {
    if (preconditioner_->RangeMap().HasComponent("face")) {
      // get upwind conductivity data
      Teuchos::RCP<CompositeVector> duw_cond = S->GetPtrW<CompositeVector>(
        Keys::getDerivKey(Keys::getKey(domain_, "upwind_overland_conductivity"), key_), name_);
      duw_cond->PutScalar(0.);

      // Then upwind.  This overwrites the boundary if upwinding says so.
      upwinding_dkdp_->Update(S);
      duw_cond->ScatterMasterToGhosted("face");
    } else {
      dcond->ScatterMasterToGhosted("cell");
    }
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void
OverlandFlow::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating BCs." << std::endl;

  AmanziMesh::Entity_ID_List cells;
  const Epetra_MultiVector& elevation =
    *S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "elevation"))->ViewComponent("face", false);

  // initialize all as null
  for (unsigned int n = 0; n != bc_markers().size(); ++n) {
    bc_markers()[n] = Operators::OPERATOR_BC_NONE;
    bc_values()[n] = 0.0;
  }

  // Head BCs are standard Dirichlet, plus elevation
  for (Functions::BoundaryFunction::Iterator bc = bc_head_->begin(); bc != bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers()[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_values()[f] = bc->second + elevation[0][f];
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc = bc_flux_->begin(); bc != bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers()[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_values()[f] = bc->second;
  }

  // zero gradient: grad h = 0 implies that q = -k grad z
  // -- cannot be done yet as rel perm update is done after this and is needed.
  // -- Instead zero gradient BCs are done in FixBCs methods.

  // Seepage face head boundary condition
  if (bc_seepage_head_->size() > 0) {
    S->GetEvaluator(Keys::getKey(domain_, "ponded_depth"))->HasFieldChanged(S.ptr(), name_);

    const CompositeVector& pd = *S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "ponded_depth"));
    const Epetra_MultiVector& h_c = *pd.ViewComponent("cell");
    const Epetra_MultiVector& elevation_c =
      *S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "elevation"))->ViewComponent("cell");

    if (pd.HasComponent("face")) {
      const Epetra_MultiVector& h_f = *pd.ViewComponent("face");
      for (Functions::BoundaryFunction::Iterator bc = bc_seepage_head_->begin();
           bc != bc_seepage_head_->end();
           ++bc) {
        int f = bc->first;

        cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        int c = cells[0];

        double hz_f = bc->second + elevation[0][f];
        double hz_c = h_c[0][c] + elevation_c[0][c];

        if (hz_f >= hz_c) {
          bc_markers()[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_values()[f] = 0.0;
        } else {
          bc_markers()[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_values()[f] = h_f[0][f] + elevation[0][f];
        }
      }

    } else {
      for (Functions::BoundaryFunction::Iterator bc = bc_seepage_head_->begin();
           bc != bc_seepage_head_->end();
           ++bc) {
        int f = bc->first;
        cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        int c = cells[0];

        double hz_f = bc->second + elevation[0][f];
        double hz_c = h_c[0][c] + elevation_c[0][c];

        if (hz_f >= hz_c) {
          bc_markers()[f] = Operators::OPERATOR_BC_NEUMANN;
          bc_values()[f] = 0.0;
        } else {
          bc_markers()[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_values()[f] = hz_f;
        }
      }
    }
  }

  // Critical depth boundary condition
  if (bc_critical_depth_->size() > 0) {
    S->GetEvaluator(Keys::getKey(domain_, "ponded_depth"))->HasFieldChanged(S.ptr(), name_);

    const Epetra_MultiVector& h_c =
      *S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "ponded_depth"))->ViewComponent("cell");
    const Epetra_MultiVector& nliq_c =
      *S->GetPtr<CompositeVector>("surface-molar_density_liquid")->ViewComponent("cell");
    double gz = -(*S->GetConstantVectorData("gravity", Tags::DEFAULT))[2];

    for (Functions::BoundaryFunction::Iterator bc = bc_critical_depth_->begin();
         bc != bc_critical_depth_->end();
         ++bc) {
      int f = bc->first;
      cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int c = cells[0];

      bc_markers()[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_values()[f] = sqrt(gz) * std::pow(h_c[0][c], 1.5) * nliq_c[0][c];
    }
  }

  // mark all remaining boundary conditions as zero flux conditions
  int nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_markers()[f] == Operators::OPERATOR_BC_NONE) {
      cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_markers()[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_values()[f] = 0.0;
      }
    }
  }
}


void
OverlandFlow::FixBCsForOperator_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the Operator." << std::endl;

  // Now we can safely calculate q = -k grad z for zero-gradient problems
  Teuchos::RCP<const CompositeVector> elev =
    S->GetPtrW<CompositeVector>(Keys::getKey(domain_, "elevation"));
  elev->ScatterMasterToGhosted();
  const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face", false);
  const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell", false);

  std::vector<WhetStone::DenseMatrix>& Aff = matrix_diff_->local_op()->matrices;

  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (Functions::BoundaryFunction::Iterator bc = bc_zero_gradient_->begin();
       bc != bc_zero_gradient_->end();
       ++bc) {
    int f = bc->first;

    AmanziMesh::Entity_ID_List cells;
    cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    AMANZI_ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];

    if (f < nfaces_owned) {
      double dp = elevation_f[0][f] - elevation_c[0][c];
      double bc_val = -dp * Aff[f](0, 0);

      bc_markers()[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_values()[f] = bc_val / mesh_->getFaceArea(f);
    }
  }
};

} // namespace Flow
} // namespace Amanzi
