/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_FECrsGraph.h"
#include "LinearOperatorFactory.hh"
#include "MatrixMFD_Coupled_Surf.hh"
#include "MatrixMFD_Surf.hh"

namespace Amanzi {
namespace Operators {

MatrixMFD_Coupled_Surf::MatrixMFD_Coupled_Surf(Teuchos::ParameterList& plist,
                                               const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
  : MatrixMFD_Coupled(plist, mesh)
{
  dump_schur_ = MatrixMFD_Coupled::dump_schur_;
  MatrixMFD_Coupled::dump_schur_ = false;
}

MatrixMFD_Coupled_Surf::MatrixMFD_Coupled_Surf(const MatrixMFD_Coupled_Surf& other)
  : MatrixMFD_Coupled(other),
    surface_mesh_(other.surface_mesh_),
    surface_A_(other.surface_A_),
    surface_B_(other.surface_B_)
{
  dump_schur_ = MatrixMFD_Coupled::dump_schur_;
  MatrixMFD_Coupled::dump_schur_ = false;
}


void
MatrixMFD_Coupled_Surf::SetSurfaceOperators(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                                            const Teuchos::RCP<MatrixMFD_TPFA>& surface_B)
{
  surface_A_ = surface_A;
  surface_B_ = surface_B;
  ASSERT(surface_A_->Mesh() == surface_B_->Mesh());
  surface_mesh_ = surface_A_->Mesh();

  MarkLocalMatricesAsChanged_();
}


void
MatrixMFD_Coupled_Surf::SetOffDiagonals(const Teuchos::RCP<const Epetra_MultiVector>& Ccc,
                                        const Teuchos::RCP<const Epetra_MultiVector>& Dcc,
                                        const Teuchos::RCP<const Epetra_MultiVector>& Ccc_surf,
                                        const Teuchos::RCP<const Epetra_MultiVector>& Dcc_surf,
                                        double scaling)
{
  MatrixMFD_Coupled::SetOffDiagonals(Ccc, Dcc, scaling);

  if (Ccc_surf == Teuchos::null) {
    Teuchos::RCP<Epetra_MultiVector> Ccc_s =
      Teuchos::rcp(new Epetra_MultiVector(surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1));
    Ccc_s->PutScalar(0.);
    Ccc_surf_ = Ccc_s;
  } else {
    Ccc_surf_ = Ccc_surf;
  }

  if (Dcc_surf == Teuchos::null) {
    Teuchos::RCP<Epetra_MultiVector> Dcc_s =
      Teuchos::rcp(new Epetra_MultiVector(surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1));
    Dcc_s->PutScalar(0.);
    Dcc_surf_ = Dcc_s;
  } else {
    Dcc_surf_ = Dcc_surf;
  }
}


void
MatrixMFD_Coupled_Surf::SymbolicAssembleGlobalMatrices()
{
  // Identical to MatrixMFD_Coupled, but does the FillMatrixGraph call
  // with a _Surf matrix, which must first be created.
  int ierr(0);
  const Epetra_BlockMap& cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL,false);
  const Epetra_BlockMap& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,false);
  const Epetra_BlockMap& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);

  // Make the double maps
  double_fmap_ = Teuchos::rcp(new Epetra_BlockMap(fmap.NumGlobalPoints(),
                                                  fmap.NumMyPoints(),
                                                  fmap.MyGlobalElements(),
                                                  2,
                                                  fmap.IndexBase(),
                                                  fmap.Comm()));

  double_cmap_ = Teuchos::rcp(new Epetra_BlockMap(cmap.NumGlobalPoints(),
                                                  cmap.NumMyPoints(),
                                                  cmap.MyGlobalElements(),
                                                  2,
                                                  cmap.IndexBase(),
                                                  cmap.Comm()));

  double_fmap_wghost_ = Teuchos::rcp(new Epetra_BlockMap(fmap_wghost.NumGlobalPoints(),
                                                         fmap_wghost.NumMyPoints(),
                                                         fmap_wghost.MyGlobalElements(),
                                                         2,
                                                         fmap_wghost.IndexBase(),
                                                         fmap_wghost.Comm()));

  // Make the matrix graphs
  int avg_entries_row = 6;
  Teuchos::RCP<Epetra_CrsGraph> cf_graph = Teuchos::rcp(
    new Epetra_CrsGraph(Copy, *double_cmap_, *double_fmap_wghost_, avg_entries_row, false));
  Teuchos::RCP<Epetra_FECrsGraph> ff_graph =
    Teuchos::rcp(new Epetra_FECrsGraph(Copy, *double_fmap_, 2 * avg_entries_row - 1, false));

  // Fill the graphs from the symbolic patterns in matrix_mfd from one of the
  // sub-block matrices.

  // IT IS ASSUMED that the two sub-block matrices fill in identical ways!
  // Create a MatrixMFD_Surf operator for filling.
  MatrixMFD_Surf blockA_surf(blockA_->plist_, blockA_->mesh_);
  blockA_surf.SetSurfaceOperator(surface_A_);
  blockA_surf.FillMatrixGraphs_(cf_graph.ptr(), ff_graph.ptr());

  // Assemble the graphs
  ierr = cf_graph->FillComplete(*double_fmap_, *double_cmap_);
  ASSERT(!ierr);
  ierr = ff_graph->GlobalAssemble(); // Symbolic graph is complete.
  ASSERT(!ierr);

  // Create the matrices
  P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));
  ierr = P2f2f_->GlobalAssemble();
  A2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));
  ierr = A2f2f_->GlobalAssemble();
  ASSERT(!ierr);
}

void
MatrixMFD_Coupled_Surf::AssembleSchur_() const
{
  // Base ComputeSchurComplement() gets the standard face parts
  MatrixMFD_Coupled::AssembleSchur_();

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& surf_cmap_wghost = surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,true);
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);

  // Grab the surface operators' TPFA cell-cell matrices
  const Epetra_FECrsMatrix& App = *surface_A_->TPFA();
  const Epetra_FECrsMatrix& Bpp = *surface_B_->TPFA();

  int entriesA = 0;
  int entriesB = 0;

  double *valuesA, *valuesB;
  valuesA = new double[9];
  valuesB = new double[9];

  int *indicesA, *indicesB;
  indicesA = new int[9];
  indicesB = new int[9];

  int subsurf_entries = 0;
  int* gsubsurfindices;
  gsubsurfindices = new int[9];
  double* gsubsurfvalues;
  gsubsurfvalues = new double[9];

  Epetra_SerialDenseMatrix block(2, 2);

  int ierr(0);

  // Now, add in the contributions from App
  // Loop over surface cells (subsurface faces)
  int ncells_surf = surface_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc = 0; sc != ncells_surf; ++sc) {
    // Access the row from the surfaces
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);

    ierr = App.ExtractGlobalRowCopy(sc_global, 9, entriesA, valuesA, indicesA);
    ASSERT(!ierr);
    ierr = Bpp.ExtractGlobalRowCopy(sc_global, 9, entriesB, valuesB, indicesB);
    ASSERT(!ierr);

    // ensure consistency of the sparsity structure, this can likely be
    // removed eventually.
    ASSERT(entriesA == entriesB);
    for (int m = 0; m != entriesA; ++m) { ASSERT(indicesA[m] == indicesB[m]); }

    // Convert local cell numbers to domain's local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
    AmanziMesh::Entity_ID frow_global = fmap_wghost.GID(frow);
    int diag = -1;
    for (int m = 0; m != entriesA; ++m) {
      indicesA[m] = surf_cmap_wghost.LID(indicesA[m]);
      if (sc == indicesA[m]) diag = m; // save the diagonal index
      indicesA[m] = surface_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, indicesA[m]);
      indicesA[m] = fmap_wghost.GID(indicesA[m]);
    }

    // Add the entries
    ierr = P2f2f_->BeginSumIntoGlobalValues(frow_global, entriesA, indicesA);
    ASSERT(!ierr);


    for (int m = 0; m != entriesA; ++m) {
      block(0, 0) = valuesA[m];
      //      block(0,0) = std::max(valuesA[m], 1.e-10);
      block(1, 1) = valuesB[m];

      if (frow_global == indicesA[m]) {
        block(0, 1) = (*Ccc_surf_)[0][sc] * scaling_;
        block(1, 0) = (*Dcc_surf_)[0][sc] * scaling_;
      } else {
        block(0, 1) = 0.;
        block(1, 0) = 0.;
      }
      ierr = P2f2f_->SubmitBlockEntry(block.A(), block.LDA(), block.M(), block.N());
      ASSERT(!ierr);
    }

    ierr = P2f2f_->EndSubmitEntries();
    ASSERT(!ierr);
  }

  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);

  // DEBUG dump
  if (dump_schur_) {
    std::stringstream filename_s;
    filename_s << "schur_MatrixMFD_Coupled_Surf_" << 0 << ".txt";
    int ierr = EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *P2f2f_);
    ASSERT(!ierr);
  }

  delete[] indicesA;
  delete[] indicesB;
  delete[] valuesA;
  delete[] valuesB;
  delete[] gsubsurfindices;
  delete[] gsubsurfvalues;
}


void
MatrixMFD_Coupled_Surf::AssembleAff_() const
{
  // Base ComputeSchurComplement() gets the standard face parts
  MatrixMFD_Coupled::AssembleAff_();

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& surf_cmap_wghost = surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,true);
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);

  // Grab the surface operators' TPFA cell-cell matrices
  const Epetra_FECrsMatrix& App = *surface_A_->TPFA();
  const Epetra_FECrsMatrix& Bpp = *surface_B_->TPFA();

  int entriesA = 0;
  int entriesB = 0;

  double *valuesA, *valuesB;
  valuesA = new double[9];
  valuesB = new double[9];

  int *indicesA, *indicesB;
  indicesA = new int[9];
  indicesB = new int[9];

  int subsurf_entries = 0;
  int* gsubsurfindices;
  gsubsurfindices = new int[9];
  double* gsubsurfvalues;
  gsubsurfvalues = new double[9];

  Epetra_SerialDenseMatrix block(2, 2);

  int ierr(0);

  // Now, add in the contributions from App
  // Loop over surface cells (subsurface faces)
  int ncells_surf = surface_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc = 0; sc != ncells_surf; ++sc) {
    // Access the row from the surfaces
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);

    ierr = App.ExtractGlobalRowCopy(sc_global, 9, entriesA, valuesA, indicesA);
    ASSERT(!ierr);
    ierr = Bpp.ExtractGlobalRowCopy(sc_global, 9, entriesB, valuesB, indicesB);
    ASSERT(!ierr);

    // ensure consistency of the sparsity structure, this can likely be
    // removed eventually.
    ASSERT(entriesA == entriesB);
    for (int m = 0; m != entriesA; ++m) { ASSERT(indicesA[m] == indicesB[m]); }

    // Convert local cell numbers to domain's local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
    AmanziMesh::Entity_ID frow_global = fmap_wghost.GID(frow);
    int diag = -1;
    for (int m = 0; m != entriesA; ++m) {
      indicesA[m] = surf_cmap_wghost.LID(indicesA[m]);
      if (sc == indicesA[m]) diag = m; // save the diagonal index
      indicesA[m] = surface_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, indicesA[m]);
      indicesA[m] = fmap_wghost.GID(indicesA[m]);
    }

    // Add the entries
    ierr = A2f2f_->BeginSumIntoGlobalValues(frow_global, entriesA, indicesA);
    ASSERT(!ierr);


    for (int m = 0; m != entriesA; ++m) {
      block(0, 0) = valuesA[m];
      //      block(0,0) = std::max(valuesA[m], 1.e-10);
      block(1, 1) = valuesB[m];

      if (frow_global == indicesA[m]) {
        block(0, 1) = (*Ccc_surf_)[0][sc] * scaling_;
        block(1, 0) = (*Dcc_surf_)[0][sc] * scaling_;
      } else {
        block(0, 1) = 0.;
        block(1, 0) = 0.;
      }
      ierr = A2f2f_->SubmitBlockEntry(block.A(), block.LDA(), block.M(), block.N());
      ASSERT(!ierr);
    }

    ierr = A2f2f_->EndSubmitEntries();
    ASSERT(!ierr);
  }

  ierr = A2f2f_->GlobalAssemble();
  ASSERT(!ierr);

  delete[] indicesA;
  delete[] indicesB;
  delete[] valuesA;
  delete[] valuesB;
  delete[] gsubsurfindices;
  delete[] gsubsurfvalues;
}


int
MatrixMFD_Coupled_Surf::Apply(const TreeVector& X, TreeVector& Y) const
{
  // Apply the coupled Matrix, which gets the subsurface.
  int ierr = MatrixMFD_Coupled::Apply(X, Y);
  ASSERT(!ierr);

  // Add in the Surface components
  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->Data();
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->Data();

  // Manually copy data -- TRILINOS FAIL
  const Epetra_MultiVector& XA_f = *XA->ViewComponent("face", false);
  const Epetra_MultiVector& XB_f = *XB->ViewComponent("face", false);
  Epetra_MultiVector surf_XA(surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
  Epetra_MultiVector surf_XB(surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1);

  for (int sc = 0; sc != surf_XA.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
    surf_XA[0][sc] = XA_f[0][f];
    surf_XB[0][sc] = XB_f[0][f];
  }

  // Apply the surface-only operators, blockwise
  Epetra_MultiVector surf_YA(surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
  Epetra_MultiVector surf_YB(surface_mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), 1);

  // -- A_surf * lambda_p ...
  ierr |= surface_A_->Apply(surf_XA, surf_YA);

  // -- + Ccc_surf * lambda_T
  ierr |= surf_YA.Multiply(scaling_, *Ccc_surf_, surf_XB, 1.0);

  // -- B_surf * lambda_T ...
  ierr |= surface_B_->Apply(surf_XB, surf_YB);

  // -- + Dcc_surf * lambda_p
  ierr |= surf_YB.Multiply(scaling_, *Dcc_surf_, surf_XA, 1.0);
  ASSERT(!ierr);

  // Add back into Y
  Epetra_MultiVector& YA_f = *YA->ViewComponent("face", false);
  Epetra_MultiVector& YB_f = *YB->ViewComponent("face", false);
  for (int sc = 0; sc != surf_YA.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
    YA_f[0][f] += surf_YA[0][sc];
    YB_f[0][f] += surf_YB[0][sc];
  }
  return ierr;
}

/* ******************************************************************
 * Solve the bottom row of the block system for lambda, given p.
 ****************************************************************** */
void
MatrixMFD_Coupled_Surf::UpdateConsistentFaceCorrection(const TreeVector& u,
                                                       const Teuchos::Ptr<TreeVector>& Pu)
{
  if (!assembled_operator_) AssembleAff_();

  // Aff solutions
  if (A2f2f_solver_ == Teuchos::null) {
    if (plist_.isSublist("consistent face solver")) {
      Teuchos::ParameterList Aff_plist = plist_.sublist("consistent face solver");
      A2f2f_op_ = Teuchos::rcp(new EpetraMatrixDefault<Epetra_FEVbrMatrix>(Aff_plist));
      A2f2f_op_->Update(A2f2f_);

      if (Aff_plist.isParameter("iterative method")) {
        AmanziSolvers::LinearOperatorFactory<EpetraMatrix, Epetra_Vector, Epetra_BlockMap> op_fac;
        A2f2f_solver_ = op_fac.Create(Aff_plist, A2f2f_op_);
      } else {
        A2f2f_solver_ = A2f2f_op_;
      }
    } else {
      Errors::Message msg("MatrixMFD::UpdateConsistentFaceConstraints was called, but no "
                          "consistent face solver sublist was provided.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // rhs = u_f - Afc * Pu_c
  Teuchos::RCP<TreeVector> work = Teuchos::rcp(new TreeVector(*Pu));
  blockA_->ApplyAfc(*Pu->SubVector(0)->Data(),
                    *work->SubVector(0)->Data(),
                    0.); // Afc is kept in the transpose form.
  blockB_->ApplyAfc(*Pu->SubVector(1)->Data(),
                    *work->SubVector(1)->Data(),
                    0.); // Afc is kept in the transpose form.
  work->SubVector(0)
    ->Data()
    ->ViewComponent("face", false)
    ->Update(1.0, *u.SubVector(0)->Data()->ViewComponent("face", false), -1.0);
  work->SubVector(1)
    ->Data()
    ->ViewComponent("face", false)
    ->Update(1.0, *u.SubVector(1)->Data()->ViewComponent("face", false), -1.0);

  // reorder rhs_f by 2xfaces
  Epetra_Vector Xf(*double_fmap_);
  Epetra_Vector Yf(*double_fmap_);
  {
    const Epetra_MultiVector& Af = *work->SubVector(0)->Data()->ViewComponent("face", false);
    const Epetra_MultiVector& Bf = *work->SubVector(1)->Data()->ViewComponent("face", false);

    for (int f = 0; f != Af.MyLength(); ++f) {
      Xf[2 * f] = Af[0][f];
      Xf[2 * f + 1] = Bf[0][f];
    }
  }

  // solve
  A2f2f_op_->Destroy();
  A2f2f_op_->Update(A2f2f_);
  int ierr = A2f2f_solver_->ApplyInverse(Xf, Yf);
  ASSERT(!ierr);

  // copy back into subblock
  {
    Epetra_MultiVector& Af = *Pu->SubVector(0)->Data()->ViewComponent("face", false);
    Epetra_MultiVector& Bf = *Pu->SubVector(1)->Data()->ViewComponent("face", false);

    for (int f = 0; f != Af.MyLength(); ++f) {
      Af[0][f] = Yf[2 * f];
      Bf[0][f] = Yf[2 * f + 1];
    }
  }
}


} // namespace Operators
} // namespace Amanzi
