/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for a StrongMPC which uses a preconditioner in which the
   block-diagonal cell-local matrix is dense.  If the system looks something
   like:

   A( y1, y2, x, t ) = 0
   B( y1, y2, x, t ) = 0

   where y1,y2 are spatially varying unknowns that are discretized using the MFD
   method (and therefore have both cell and face unknowns), an approximation to
   the Jacobian is written as

   [  dA_c/dy1_c  dA_c/dy1_f   dA_c/dy2_c       0      ]
   [  dA_f/dy1_c  dA_f/dy1_f      0              0      ]
   [  dB_c/dy1_c     0          dB_c/dy2_c  dB_c/dy2_f ]
   [      0           0          dB_f/dy2_c  dB_f/dy2_f ]


   Note that the upper left block is the standard preconditioner for the A
   system, and the lower right block is the standard precon for the B system,
   and we have simply added cell-based couplings, dA_c/dy2_c and dB_c/dy1_c.

   In the temperature/pressure system, these correspond to d_water_content /
   d_temperature and d_energy / d_pressure.

   ------------------------------------------------------------------------- */

#include <fstream>
#include "EpetraExt_RowMatrixOut.h"

#include "FieldEvaluator.hh"
#include "Operator.hh"
#include "TreeOperator.hh"
#include "PDE_Accumulation.hh"

#include "mpc_coupled_cells.hh"

namespace Amanzi {

void MPCCoupledCells::Setup(const Teuchos::Ptr<State>& S) {
  StrongMPC<PK_PhysicalBDF_Default>::Setup(S);

  A_key_ = plist_->get<std::string>("conserved quantity A");
  B_key_ = plist_->get<std::string>("conserved quantity B");
  y1_key_ = plist_->get<std::string>("primary variable A");
  y2_key_ = plist_->get<std::string>("primary variable B");
  dA_dy2_key_ = std::string("d")+A_key_+std::string("_d")+y2_key_;
  dB_dy1_key_ = std::string("d")+B_key_+std::string("_d")+y1_key_;

  Key mesh_key = plist_->get<std::string>("domain name");
  mesh_ = S->GetMesh(mesh_key);

  // set up debugger
  db_ = sub_pks_[0]->debugger();

  // Get the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Operator> pcA = sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Operator> pcB = sub_pks_[1]->preconditioner();

  // Create the combined operator
  Teuchos::RCP<TreeVectorSpace> tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcA->DomainMap()))));
  tvs->PushBack(Teuchos::rcp(new TreeVectorSpace(Teuchos::rcpFromRef(pcB->DomainMap()))));

  preconditioner_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  preconditioner_->set_operator_block(0, 0, pcA);
  preconditioner_->set_operator_block(1, 1, pcB);
  
  // create coupling blocks and push them into the preconditioner...
  S->RequireField(A_key_)->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(y2_key_)->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(A_key_);
  S->RequireFieldEvaluator(y2_key_);

  if (!plist_->get<bool>("no dA/dy2 block", false)) {
    Teuchos::ParameterList& acc_pc_plist = plist_->sublist("dA_dy2 accumulation preconditioner");
    acc_pc_plist.set("entity kind", "cell");
    dA_dy2_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, mesh_));
    preconditioner_->set_operator_block(0, 1, dA_dy2_->global_operator());
  }

  S->RequireField(B_key_)->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(y1_key_)->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(B_key_);
  S->RequireFieldEvaluator(y1_key_);

  if (!plist_->get<bool>("no dB/dy1 block", false)) {
    Teuchos::ParameterList& acc_pc_plist = plist_->sublist("dB_dy1 accumulation preconditioner");
    acc_pc_plist.set("entity kind", "cell");
    dB_dy1_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, mesh_));
    preconditioner_->set_operator_block(1, 0, dB_dy1_->global_operator());
  }

  // setup and initialize the preconditioner
  preconditioner_->set_inverse_parameters(plist_->sublist("preconditioner"));
}


// updates the preconditioner
void MPCCoupledCells::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
  StrongMPC<PK_PhysicalBDF_Default>::UpdatePreconditioner(t,up,h);

  if (dA_dy2_ != Teuchos::null &&
      S_next_->GetFieldEvaluator(A_key_)->IsDependency(S_next_.ptr(), y2_key_)) {
    dA_dy2_->global_operator()->Init();
    S_next_->GetFieldEvaluator(A_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y2_key_);
    Teuchos::RCP<const CompositeVector> dA_dy2_v = S_next_->GetFieldData(dA_dy2_key_);
    db_->WriteVector("  dwc_dT", dA_dy2_v.ptr());
    
    // -- update the cell-cell block
    dA_dy2_->AddAccumulationTerm(*dA_dy2_v, h, "cell", false);
  }

  if (dB_dy1_ != Teuchos::null &&
      S_next_->GetFieldEvaluator(B_key_)->IsDependency(S_next_.ptr(), y1_key_)) {
    dB_dy1_->global_operator()->Init();
    S_next_->GetFieldEvaluator(B_key_)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y1_key_);
    Teuchos::RCP<const CompositeVector> dB_dy1_v = S_next_->GetFieldData(dB_dy1_key_);
    db_->WriteVector("  dE_dp", dB_dy1_v.ptr());
    
    // -- update the cell-cell block
    dB_dy1_->AddAccumulationTerm(*dB_dy1_v, h, "cell", false);
  }

}


// applies preconditioner to u and returns the result in Pu
int MPCCoupledCells::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Residuals:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  r_p"); vnames.push_back("  r_T"); 
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(u->SubVector(0)->Data().ptr()); 
    vecs.push_back(u->SubVector(1)->Data().ptr()); 
    db_->WriteVectors(vnames, vecs, true);
  }
  
  int ierr = preconditioner_->ApplyInverse(*u, *Pu);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "PC * residuals:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  PC*r_p"); vnames.push_back("  PC*r_T"); 
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(Pu->SubVector(0)->Data().ptr()); 
    vecs.push_back(Pu->SubVector(1)->Data().ptr()); 
    db_->WriteVectors(vnames, vecs, true);
  }

  return (ierr > 0) ? 0 : 1;
}


} //  namespace
