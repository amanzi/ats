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
void PK_PhysicalBDF_Default::Setup(const Teuchos::Ptr<State>& S) {

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
void PK_PhysicalBDF_Default::Initialize(const Teuchos::Ptr<State>& S) {
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
        Teuchos::RCP<const TreeVector> res) {
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

    if (*comp == std::string("cell")) {
      // error done relative to extensive, conserved quantity
      int ncells = dvec->size(*comp,false);
      for (unsigned int c=0; c!=ncells; ++c) {
        double enorm_c = std::abs(h * dvec_v[0][c])
            / (atol_*cv[0][c] + rtol_*std::abs(conserved[0][c]));

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


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void
PK_PhysicalBDF_Default::ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u) {
  auto& markers = bc_markers();
  auto& values = bc_values();
  if (u->HasComponent("face")) {
    Epetra_MultiVector& u_f = *u->ViewComponent("face",false);
    unsigned int nfaces = u_f.MyLength();
    for (unsigned int f=0; f!=nfaces; ++f) {
      if (markers[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_f[0][f] = values[f];
      }
    }
  } else if (u->HasComponent("boundary_face")) {
    const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
    const Epetra_Map& face_map = mesh_->face_map(false);

    Epetra_MultiVector& u_bf = *u->ViewComponent("boundary_face",false);
    unsigned int nfaces = u_bf.MyLength();
    for (unsigned int bf=0; bf!=nfaces; ++bf) {
      AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
      if (markers[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_bf[0][bf] = values[f];
      }
    }
  }    
};


double PK_PhysicalBDF_Default::BoundaryValue(const Teuchos::RCP<const Amanzi::CompositeVector>& solution, int face_id){
  double value=0.;

  // if (solution->HasComponent("face")){
  //   const Epetra_MultiVector& u = *solution -> ViewComponent("face",false);
  //   value = u[0][face_id];
  // }
  // else if  (solution->HasComponent("boundary_face")){
  //   const Epetra_MultiVector& u = *solution -> ViewComponent("boundary_face",false);
  //   const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  //   const Epetra_Map& f_map = mesh_->face_map(false);

  //   int face_gid = f_map.GID(face_id);
  //   int face_lbid = fb_map.LID(face_gid);

  //   value =  u[0][face_lbid];
  // }
  // else{
  //   Errors::Message msg("No component is defined for boundary faces\n");
  //   Exceptions::amanzi_throw(msg);
  // }

  // return value;
  if (solution->HasComponent("face")){
    const Epetra_MultiVector& u = *solution->ViewComponent("face",false);
    value = u[0][face_id];
  // } else if  (solution->HasComponent("boundary_face") &&
  //             bc_markers_[face_id] == Operators::OPERATOR_BC_DIRICHLET){
  } else if (bc_markers()[face_id] == Operators::OPERATOR_BC_DIRICHLET) {
    // const Epetra_MultiVector& u = *solution->ViewComponent("boundary_face",false);
    // const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
    // const Epetra_Map& f_map = mesh_->face_map(false);

    // int face_gid = f_map.GID(face_id);
    // int face_lbid = fb_map.LID(face_gid);

    // value = u[0][face_lbid];
    value = bc_values()[face_id];
  } else {
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(face_id, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
    const Epetra_MultiVector& u = *solution->ViewComponent("cell",false);
    value = u[0][cells[0]];    
  }
  return value;
}


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


int PK_PhysicalBDF_Default::BoundaryDirection(int face_id) {
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(face_id, AmanziMesh::Parallel_type::ALL, &cells);
  AMANZI_ASSERT(cells.size() == 1);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(cells[0], &faces, &dirs);
  return dirs[std::find(faces.begin(), faces.end(), face_id) - faces.begin()];
}

// Experimental approach -- calling this indicates that the time
// integration scheme is changing the value of the solution in
// state.
// -----------------------------------------------------------------------------
void PK_PhysicalBDF_Default::ChangedSolution(const Teuchos::Ptr<State>& S) {
  if (S == Teuchos::null) {
    solution_evaluator_->SetFieldAsChanged(S_next_.ptr());

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
// the value of the solution in state.
// -----------------------------------------------------------------------------
void PK_PhysicalBDF_Default::ChangedSolution() {
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
};
  

} // namespace
