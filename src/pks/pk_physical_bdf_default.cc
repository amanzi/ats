/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Standard base for most PKs, this combines both domains/meshes of
PKPhysicalBase and BDF methods of PK_BDF_Default.
------------------------------------------------------------------------- */

#include "boost/math/special_functions/fpclassify.hpp"

#include "pk_physical_bdf_default.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void PK_PhysicalBDF_Default::Setup(const Teuchos::Ptr<State>& S)
{
  // call the meat of the base constructurs via Setup methods
  PK_Physical_Default::Setup(S);
  PK_BDF_Default::Setup(S);

  // boundary conditions
  bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // convergence criteria is based on a conserved quantity
  if (conserved_key_.empty()) {
    conserved_key_ = Keys::readKey(*plist_, domain_, "conserved quantity");
  }
  S->RequireField(conserved_key_)->SetMesh(mesh_)
      ->AddComponent("cell",AmanziMesh::CELL,true);
  S->RequireFieldEvaluator(conserved_key_);

  // cell volume used throughout
  if (cell_vol_key_.empty()) {
    cell_vol_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
  }
  S->RequireField(cell_vol_key_)->SetMesh(mesh_)
      ->AddComponent("cell",AmanziMesh::CELL,true);
  S->RequireFieldEvaluator(cell_vol_key_);

  atol_ = plist_->get<double>("absolute error tolerance",1.0);
  rtol_ = plist_->get<double>("relative error tolerance",1.0);
  fluxtol_ = plist_->get<double>("flux error tolerance",1.0);
};


// -----------------------------------------------------------------------------
// initialize.  Note both BDFBase and PhysicalBase have initialize()
// methods, so we need a unique overrider.
// -----------------------------------------------------------------------------
void PK_PhysicalBDF_Default::Initialize(const Teuchos::Ptr<State>& S)
{
  // Just calls both subclass's initialize.  NOTE - order is important here --
  // PhysicalBase grabs the primary variable and stuffs it into the solution,
  // which must be done prior to BDFBase initializing the timestepper.
  PK_Physical_Default::Initialize(S);
  PK_BDF_Default::Initialize(S);
}


// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double PK_PhysicalBDF_Default::ErrorNorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> res)
{
  // Abs tol based on old conserved quantity -- we know these have been vetted
  // at some level whereas the new quantity is some iterate, and may be
  // anything from negative to overflow.
  S_inter_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& conserved = *S_inter_->GetFieldData(conserved_key_)
      ->ViewComponent("cell",true);
  const Epetra_MultiVector& cv = *S_inter_->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",true);

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << conserved_key_ << ": " << std::endl;

  Teuchos::RCP<const CompositeVector> dvec = res->Data();
  double h = S_next_->time() - S_inter_->time();

  Teuchos::RCP<const Comm_type> comm_p = mesh_->get_comm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp=dvec->begin();
       comp!=dvec->end(); ++comp) {
    double enorm_comp = 0.0;
    int enorm_loc = -1;
    const Epetra_MultiVector& dvec_v = *dvec->ViewComponent(*comp, false);

    if (*comp == "cell") {
      // error done relative to extensive, conserved quantity
      int ncells = dvec->size(*comp,false);
      for (unsigned int c=0; c!=ncells; ++c) {
        double enorm_c = std::abs(h * dvec_v[0][c])
            / (atol_*cv[0][c] + rtol_*std::abs(conserved[0][c]));
        AMANZI_ASSERT((atol_*cv[0][c] + rtol_*std::abs(conserved[0][c])) > 0.);

        if (enorm_c > enorm_comp) {
          enorm_comp = enorm_c;
          enorm_loc = c;
        }
      }

    } else if (*comp == std::string("face")) {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec->size(*comp, false);

      for (unsigned int f=0; f!=nfaces; ++f) {
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
        double cv_min = cells.size() == 1 ? cv[0][cells[0]]
            : std::min(cv[0][cells[0]],cv[0][cells[1]]);
        double conserved_min = cells.size() == 1 ? conserved[0][cells[0]]
            : std::min(conserved[0][cells[0]],conserved[0][cells[1]]);

        double enorm_f = fluxtol_ * h * std::abs(dvec_v[0][f])
            / (atol_*cv_min + rtol_*std::abs(conserved_min));
        AMANZI_ASSERT((atol_*cv_min + rtol_*std::abs(conserved_min)) > 0.);
        if (enorm_f > enorm_comp) {
          enorm_comp = enorm_f;
          enorm_loc = f;
        }
      }

    } else {
      // double norm;
      // dvec_v.Norm2(&norm);
      //      AMANZI_ASSERT(norm < 1.e-15);
    }

    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm(0.);
      dvec_v.NormInf(&infnorm);

      ENorm_t err;
      ENorm_t l_err;
      l_err.value = enorm_comp;
      l_err.gid = dvec_v.Map().GID(enorm_loc);

      int ierr;

      ierr = MPI_Allreduce(&l_err, &err, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
      AMANZI_ASSERT(!ierr);
      *vo_->os() << "  ENorm (" << *comp << ") = " << err.value << "[" << err.gid << "] (" << infnorm << ")" << std::endl;
    }

    enorm_val = std::max(enorm_val, enorm_comp);
  }

  double enorm_val_l = enorm_val;

  int ierr;
  ierr = MPI_Allreduce(&enorm_val_l, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, comm);
  AMANZI_ASSERT(!ierr);
  return enorm_val;
};


  // void PK_PhysicalBDF_Default::Solution_to_State(TreeVector& solution,
  //                                                 const Teuchos::RCP<State>& S){
  //   PK_Physical_Default::Solution_to_State(solution, S);
  // }

  // void PK_PhysicalBDF_Default::Solution_to_State(const TreeVector& soln,
  //                                                const Teuchos::RCP<State>& S){
  //   TreeVector* soln_nc_ptr = const_cast<TreeVector*>(&soln);
  //   PK_Physical_Default::Solution_to_State(soln_nc_ptr, S);
  // }

void PK_PhysicalBDF_Default::set_states(const Teuchos::RCP<State>& S,
                                        const Teuchos::RCP<State>& S_inter,
                                        const Teuchos::RCP<State>& S_next) {
  PK_Physical_Default::set_states(S, S_inter, S_next);
}


// -----------------------------------------------------------------------------
// Calling this indicates that the time integration scheme is changing the
// value of the solution in state.
// -----------------------------------------------------------------------------
void PK_PhysicalBDF_Default::ChangedSolution(const Teuchos::Ptr<State>& S) {
  if (S == Teuchos::null) {
    solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
  } else if (S == S_next_.ptr()) {
    if (!solution_evaluator_.get()) {
      Teuchos::RCP<FieldEvaluator> fm = S_next_->GetFieldEvaluator(key_);
      solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
      AMANZI_ASSERT(solution_evaluator_ != Teuchos::null);
    }
    solution_evaluator_->SetFieldAsChanged(S);
  } else {
    Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator(key_);
    Teuchos::RCP<PrimaryVariableFieldEvaluator> solution_evaluator =
      Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
    AMANZI_ASSERT(solution_evaluator != Teuchos::null);
    solution_evaluator->SetFieldAsChanged(S);
  }
};


// -----------------------------------------------------------------------------
// Calling this indicates that the time integration scheme is changing
// the value of the solution in state.  This is the variant called by the TI.
// -----------------------------------------------------------------------------
void PK_PhysicalBDF_Default::ChangedSolution() {
  solution_evaluator_->SetFieldAsChanged(Teuchos::null);
};

} // namespace
