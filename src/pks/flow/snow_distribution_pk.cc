/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"
#include "Point.hh"

#include "FunctionFactory.hh"
#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"
#include "EvaluatorPrimary.hh"

#include "PDE_DiffusionFV.hh"
#include "upwind_potential_difference.hh"
#include "upwind_total_flux.hh"

#include "snow_distribution.hh"

#include "pk_helpers.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
SnowDistribution::SnowDistribution(Teuchos::ParameterList& pk_tree,
                                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, plist, S, solution),
    PK_PhysicalBDF_Default(pk_tree, plist, S, solution),
    my_next_time_(-9.e80)
{
  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .01); // h * nl

  // Keys
  cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
  cond_key_ = Keys::readKey(*plist_, domain_, "conductivity", "conductivity");
  elev_key_ = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  precip_key_ = Keys::readKey(*plist_, domain_, "precipitation", "precipitation");
  uw_cond_key_ = Keys::readKey(*plist_, domain_, "upwind conductivity", "upwind_conductivity");
  flux_dir_key_ = Keys::readKey(*plist_, domain_, "flux direction", "flux_direction");
  potential_key_ = Keys::readKey(*plist_, domain_, "skin potential", "skin_potential");
  precip_func_key_ =
    Keys::readKey(*plist_, domain_, "precipitation function", "precipitation_function");

  dt_factor_ = plist_->get<double>("distribution time", 86400.0);

  // -- elevation evaluator
  bool standalone_elev = S->GetMesh() == S->GetMesh(domain_);
  if (!standalone_elev && !S->FEList().isSublist(elev_key_)) {
    S->GetEvaluatorList(elev_key_).set("evaluator type", "meshed elevation");
  }
}


void
SnowDistribution::Setup()
{
  PK_PhysicalBDF_Default::Setup();
  SetupSnowDistribution_();
  SetupPhysicalEvaluators_();
}


void
SnowDistribution::SetupSnowDistribution_()
{
  // precip function
  Teuchos::ParameterList& precip_func = plist_->sublist("precipitation function");
  FunctionFactory fac;
  precip_func_ = Teuchos::rcp(fac.Create(precip_func));

  // -- get conserved variable (snow-precip) and evaluator and derivative for PC
  requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(conserved_key_, tag_current_, *S_, name_);

  // -- cell volume and evaluator
  S_->Require<CompositeVector, CompositeVectorSpace>(cv_key_, tag_next_)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(cv_key_, tag_next_);

  // elevation
  S_->Require<CompositeVector, CompositeVectorSpace>(elev_key_, tag_next_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  S_->RequireEvaluator(elev_key_, tag_next_);

  // boundary conditions
  auto& markers = bc_markers();
  auto& values = bc_values();
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  markers.resize(nfaces, Operators::OPERATOR_BC_NONE);
  values.resize(nfaces, 0.0);
  UpdateBoundaryConditions_(tag_next_); // never change

  // -- create the upwinding method.
  S_->Require<CompositeVector, CompositeVectorSpace>(uw_cond_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::FACE, 1);
  S_->GetRecordW(uw_cond_key_, tag_next_, name_).set_io_vis(false);

  // flux direction required for upwinding
  S_->Require<CompositeVector, CompositeVectorSpace>(flux_dir_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->SetComponent("face", AmanziMesh::FACE, 1);

  upwind_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
  upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_, tag_next_, flux_dir_key_, 1.e-8));

  // -- operator for the diffusion terms
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", "upwind: face");
  if (mfd_plist.isParameter("Newton correction")) {
    Errors::Message message;
    message << name_ << ": The forward operator for Diffusion should not set a "
            << "\"Newton correction\" term, perhaps you meant to put this in a "
            << "\"Diffusion PC\" sublist.";
    Exceptions::amanzi_throw(message);
  }
  matrix_diff_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(mfd_plist, mesh_));
  matrix_diff_->SetBCs(bc_, bc_);
  matrix_diff_->SetTensorCoefficient(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // -- create the operator, data for flux directions for upwinding
  Teuchos::ParameterList& face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(face_diff_list, mesh_));
  face_matrix_diff_->SetTensorCoefficient(Teuchos::null);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  face_matrix_diff_->SetBCs(bc_, bc_);
  face_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // -- preconditioner
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("inverse", plist_->sublist("inverse"));
  auto& inv_list = mfd_pc_plist.sublist("inverse");
  inv_list.setParameters(plist_->sublist("preconditioner"));
  inv_list.setParameters(plist_->sublist("linear solver"));

  preconditioner_diff_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(mfd_pc_plist, mesh_));
  preconditioner_diff_->SetBCs(bc_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  // accumulation operator for the preconditioenr
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ =
    Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));
}


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void
SnowDistribution::SetupPhysicalEvaluators_()
{
  // -- evaluator for potential field, h + z
  requireAtNext(potential_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- snow_conductivity evaluator
  requireAtNext(cond_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
      cond_key_, tag_next_, key_, tag_next_)
    .SetGhosted();
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void
SnowDistribution::Initialize()
{
  // Initialize BDF stuff and physical domain stuff.
  PK_PhysicalBDF_Default::Initialize();

  // Set extra fields as initialized -- these don't currently have evaluators.
  S_->GetW<CompositeVector>(uw_cond_key_, tag_next_, name_).PutScalar(1.0);
  S_->GetRecordW(uw_cond_key_, tag_next_, name_).set_initialized();

  if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    S_->GetW<CompositeVector>(flux_dir_key_, tag_next_, name_).PutScalar(0.);
    S_->GetRecordW(flux_dir_key_, tag_next_, name_).set_initialized();
  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool
SnowDistribution::UpdatePermeabilityData_(const Tag& tag)
{
  bool update_perm = S_->GetEvaluator(cond_key_, tag).Update(*S_, name_);
  update_perm |= S_->GetEvaluator(precip_key_, tag).Update(*S_, name_);
  update_perm |= S_->GetEvaluator(potential_key_, tag).Update(*S_, name_);

  if (update_perm) {
    if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<CompositeVector> flux_dir =
        S_->GetPtrW<CompositeVector>(flux_dir_key_, tag, name_);

      // Derive the flux
      Teuchos::RCP<const CompositeVector> potential =
        S_->GetPtr<CompositeVector>(potential_key_, tag);
      face_matrix_diff_->UpdateFlux(potential.ptr(), flux_dir.ptr());
    }

    // get snow_conductivity data
    Teuchos::RCP<const CompositeVector> cond = S_->GetPtr<CompositeVector>(cond_key_, tag);

    // get upwind snow_conductivity data
    Teuchos::RCP<CompositeVector> uw_cond = S_->GetPtrW<CompositeVector>(uw_cond_key_, tag, name_);

    { // place interior cells on boundary faces
      const auto& cond_c = *cond->ViewComponent("cell", false);
      auto& uw_cond_f = *uw_cond->ViewComponent("face", false);
      int nfaces = uw_cond_f.MyLength();
      AmanziMesh::Entity_ID_List cells;
      for (int f = 0; f != nfaces; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        if (cells.size() == 1) {
          int c = cells[0];
          uw_cond_f[0][f] = cond_c[0][c];
        }
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(*cond, *uw_cond, *S_);
    uw_cond->ScatterMasterToGhosted("face");
  }
  return update_perm;
}

void
SnowDistribution::UpdateBoundaryConditions_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating BCs." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  // mark all remaining boundary conditions as zero flux conditions
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List cells;
  for (int f = 0; f < nfaces_owned; f++) {
    if (markers[f] == Operators::OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
      }
    }
  }
}

} // namespace Flow
} // namespace Amanzi
