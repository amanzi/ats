/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Author: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"
#include "Point.hh"

#include "FunctionFactory.hh"
#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"
#include "independent_variable_field_evaluator.hh"

#include "PDE_DiffusionFV.hh"
#include "upwind_potential_difference.hh"
#include "upwind_total_flux.hh"

#include "snow_distribution.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1

SnowDistribution::SnowDistribution(Teuchos::ParameterList& FElist,
                                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    PK_PhysicalBDF_Default(FElist, plist, S, solution),
    my_next_time_(-9.e80)
{
  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .01); // h * nl

  dt_factor_ = plist_->get<double>("distribution time", 86400.0);
}


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void SnowDistribution::Setup(const Teuchos::Ptr<State>& S) {
  //PKPhysicalBDFBase::setup(S);
  PK_PhysicalBDF_Default::Setup(S);
  SetupSnowDistribution_(S);
  SetupPhysicalEvaluators_(S);
}


void SnowDistribution::SetupSnowDistribution_(const Teuchos::Ptr<State>& S) {
  // precip function
  Teuchos::ParameterList& precip_func = plist_->sublist("precipitation function");
  FunctionFactory fac;
  precip_func_ = Teuchos::rcp(fac.Create(precip_func));

  // Require fields and evaluators for those fields.
  S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  // -- cell volume and evaluator
  S->RequireFieldEvaluator(Keys::getKey(domain_, "cell_volume"));
  
  // boundary conditions
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  bc_markers().resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values().resize(nfaces, 0.0);
  UpdateBoundaryConditions_(S); // never change

  // -- create the upwinding method.
  S->RequireField(Keys::getKey(domain_,"upwind_conductivity"), name_)->SetMesh(mesh_)
    ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->GetField(Keys::getKey(domain_,"upwind_conductivity"),name_)->set_io_vis(false);

  upwind_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
  upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
          Keys::getKey(domain_,"conductivity"),
          Keys::getKey(domain_,"upwind_conductivity"),
          Keys::getKey(domain_,"flux_direction"), 1.e-8));

  // -- operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_->sublist("diffusion");
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
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(face_diff_list, mesh_));
  face_matrix_diff_->SetTensorCoefficient(Teuchos::null);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  face_matrix_diff_->SetBCs(bc_, bc_);
  face_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  S->RequireField(Keys::getKey(domain_, "flux_direction"), name_)
      ->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);

  // -- preconditioner
  Teuchos::ParameterList mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("inverse", plist_->sublist("inverse"));
  // old style... deprecate me!
  mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("preconditioner"));
  mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("linear solver"));
  
  preconditioner_diff_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(mfd_pc_plist, mesh_));
  preconditioner_diff_->SetBCs(bc_, bc_);
  preconditioner_diff_->SetTensorCoefficient(Teuchos::null);

  preconditioner_ = preconditioner_diff_->global_operator();
  
  // accumulation operator for the preconditioenr
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  // symbolic assemble, get PC
}


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void SnowDistribution::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- evaluator for potential field, h + z

  S->RequireFieldEvaluator(Keys::getKey(domain_,"skin_potential"));
  S->RequireField(Keys::getKey(domain_,"skin_potential"))->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- snow_conductivity evaluator
  S->RequireFieldEvaluator(Keys::getKey(domain_,"conductivity"));
  S->RequireField(Keys::getKey(domain_,"conductivity"))->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);

}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void SnowDistribution::Initialize(const Teuchos::Ptr<State>& S) {
  // Initialize BDF stuff and physical domain stuff.
  PK_PhysicalBDF_Default::Initialize(S);

  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData(Keys::getKey(domain_,"upwind_conductivity"),name_)->PutScalar(1.0);
  S->GetField(Keys::getKey(domain_,"upwind_conductivity"),name_)->set_initialized();

  if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    S->GetFieldData(Keys::getKey(domain_,"flux_direction"), name_)->PutScalar(0.);
    S->GetField(Keys::getKey(domain_,"flux_direction"), name_)->set_initialized();
  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool SnowDistribution::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {

  bool update_perm = S->GetFieldEvaluator(Keys::getKey(domain_,"conductivity"))
      ->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator(Keys::getKey(domain_,"precipitation"))->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator(Keys::getKey(domain_,"skin_potential"))->HasFieldChanged(S, name_);

  if (update_perm) {
    if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<CompositeVector> flux_dir =
          S->GetFieldData(Keys::getKey(domain_,"flux_direction"), name_);

      // Derive the flux
      Teuchos::RCP<const CompositeVector> potential = S->GetFieldData(Keys::getKey(domain_,"skin_potential"));
      face_matrix_diff_->UpdateFlux(potential.ptr(), flux_dir.ptr());
    }

    // get snow_conductivity data
    const Epetra_MultiVector& cond_c = *S->GetFieldData(Keys::getKey(domain_,"conductivity"))
        ->ViewComponent("cell",false);

    // get upwind snow_conductivity data
    Teuchos::RCP<CompositeVector> uw_cond =
      S->GetFieldData(Keys::getKey(domain_,"upwind_conductivity"), name_);

    { // place interior cells on boundary faces
      Epetra_MultiVector& uw_cond_f = *uw_cond->ViewComponent("face",false);

      int nfaces = uw_cond_f.MyLength();
      AmanziMesh::Entity_ID_List cells;
      for (int f=0; f!=nfaces; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        if (cells.size() == 1) {
          int c = cells[0];
          uw_cond_f[0][f] = cond_c[0][c];
        }
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);
    uw_cond->ScatterMasterToGhosted("face");
  }

  return update_perm;
}

void SnowDistribution::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

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

} // namespace
} // namespace

