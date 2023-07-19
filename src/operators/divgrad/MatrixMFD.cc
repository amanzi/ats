/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD.hh"
#include "LinearOperatorFactory.hh"

namespace Amanzi {
namespace Operators {

#define APPLY_UNASSEMBLED 1

/* ******************************************************************
 * Constructor
 ****************************************************************** */
MatrixMFD::MatrixMFD(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : plist_(plist),
    mesh_(mesh),
    flag_symmetry_(false),
    assembled_operator_(false),
    assembled_schur_(false),
    assembled_rhs_(false),
    method_(MFD3D_NULL)
{
  InitializeFromPList_();
}


/* ******************************************************************
 * Copy constructor
 ****************************************************************** */
MatrixMFD::MatrixMFD(const MatrixMFD& other)
  : plist_(other.plist_),
    mesh_(other.mesh_),
    flag_symmetry_(other.flag_symmetry_),
    assembled_operator_(false),
    assembled_schur_(false),
    assembled_rhs_(false),
    method_(other.method_)
{
  InitializeFromPList_();
}


/* ******************************************************************
 * operator= copies local matrices
 ****************************************************************** */
MatrixMFD&
MatrixMFD::operator=(const MatrixMFD& other)
{
  if (this != &other) {
    Mff_cells_ = other.Mff_cells_;
    Aff_cells_ = other.Aff_cells_;
    Acf_cells_ = other.Acf_cells_;
    Afc_cells_ = other.Afc_cells_;
    Ff_cells_ = other.Ff_cells_;
    Fc_cells_ = other.Fc_cells_;
  }
  return *this;
}


/* ******************************************************************
 * Initialization of method, solver, etc.
 ****************************************************************** */
void
MatrixMFD::InitializeFromPList_()
{
  std::string methodstring = plist_.get<std::string>("MFD method");
  method_ = MFD3D_NULL;

  // standard MFD
  if (methodstring == "monotone mfd hex") { // two monotone methods
    method_ = MFD3D_HEXAHEDRA_MONOTONE;
  } else if (methodstring == "monotone mfd") {
    method_ = MFD3D_POLYHEDRA_MONOTONE;
  } else if (methodstring == "support operator") {
    method_ = MFD3D_SUPPORT_OPERATOR;
  } else if (methodstring == "two point flux approximation") {
    method_ = MFD3D_TPFA;
  } else if (methodstring == "finite volume") {
    method_ = FV_TPFA;
  } else if (methodstring == "optimized mfd") {
    method_ = MFD3D_OPTIMIZED;
  } else if (methodstring == "optimized mfd scaled") {
    method_ = MFD3D_OPTIMIZED_SCALED;
  } else if (methodstring == "mfd") { // first basic mfd
    method_ = MFD3D_POLYHEDRA;
  } else if (methodstring == "mfd scaled") { // second basic mfd
    method_ = MFD3D_POLYHEDRA_SCALED;
  } else {
    Errors::Message msg("MatrixMFD: unexpected discretization method");
    Exceptions::amanzi_throw(msg);
  }

  // vector space
  std::vector<std::string> names;
  std::vector<AmanziMesh::Entity_kind> locations;
  if (method_ != FV_TPFA) {
    names.push_back("cell");
    names.push_back("face");
    locations.push_back(AmanziMesh::Entity_kind::CELL);
    locations.push_back(AmanziMesh::Entity_kind::FACE);
  } else {
    names.push_back("cell");
    names.push_back("boundary_face");
    locations.push_back(AmanziMesh::Entity_kind::CELL);
    locations.push_back(AmanziMesh::Entity_kind::BOUNDARY_FACE);
  }
  std::vector<int> ndofs(2, 1);

  space_ = Teuchos::rcp(new CompositeVectorSpace());
  space_->SetMesh(mesh_)->SetGhosted()->SetComponents(names, locations, ndofs);

  // preconditioner
  if (plist_.isSublist("preconditioner")) {
    Teuchos::ParameterList pc_list = plist_.sublist("preconditioner");
    AmanziPreconditioners::PreconditionerFactory pc_fac;
    S_pc_ = pc_fac.Create(pc_list);
    AmanziPreconditioners::PreconditionerFactory pc_fac2;
    Aff_pc_ = pc_fac2.Create(pc_list);
  }

  // verbose object
  vo_ = Teuchos::rcp(new VerboseObject("MatrixMFD", plist_));
}


// main computational methods
/* ******************************************************************
 * Calculate elemental inverse mass matrices.
 ****************************************************************** */
void
MatrixMFD::CreateMFDmassMatrices(const Teuchos::Ptr<std::vector<WhetStone::Tensor>>& K)
{
  // tag global matrices as invalid
  MarkLocalMatricesAsChanged_();

  int dim = mesh_->getSpaceDimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  //  AmanziMesh::Entity_ID_List faces;

  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  if (Mff_cells_.size() != ncells) { Mff_cells_.resize(static_cast<size_t>(ncells)); }

  int ok;
  nokay_ = npassed_ = 0;

  WhetStone::Tensor Kc;
  if (K == Teuchos::null) {
    Kc.Init(mesh_->getSpaceDimension(), 1);
    Kc(0, 0) = 1.0;
  }

  //int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);
  for (int c = 0; c != ncells; ++c) {
    int nfaces = mesh_->getCellNumFaces(c);

    WhetStone::DenseMatrix Mff(nfaces, nfaces);

    if (K != Teuchos::null) { Kc = (*K)[c]; }

    if (method_ == MFD3D_POLYHEDRA_SCALED) {
      ok = mfd.MassMatrixInverseScaled(c, Kc, Mff);
    } else if (method_ == MFD3D_POLYHEDRA_MONOTONE) {
      ok = mfd.MassMatrixInverseMMatrix(c, Kc, Mff);
      if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_WRONG) {
        ok = mfd.MassMatrixInverseTPFA(c, Kc, Mff);
        nokay_--;
        npassed_++;
      }
    } else if (method_ == MFD3D_POLYHEDRA) {
      ok = mfd.MassMatrixInverse(c, Kc, Mff);
    } else if (method_ == MFD3D_OPTIMIZED_SCALED) {
      ok = mfd.MassMatrixInverseOptimizedScaled(c, Kc, Mff);
    } else if (method_ == MFD3D_OPTIMIZED) {
      ok = mfd.MassMatrixInverseOptimized(c, Kc, Mff);
    } else if (method_ == MFD3D_HEXAHEDRA_MONOTONE) {
      if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2))
        ok = mfd.MassMatrixInverseMMatrixHex(c, Kc, Mff);
      else
        ok = mfd.MassMatrixInverse(c, Kc, Mff);
    } else if (method_ == MFD3D_TPFA) {
      ok = mfd.MassMatrixInverseTPFA(c, Kc, Mff);
    } else if (method_ == MFD3D_SUPPORT_OPERATOR) {
      ok = mfd.MassMatrixInverseSO(c, Kc, Mff);
    } else {
      Errors::Message msg(
        "MatrixMFD: unexpected discretization methods (contact lipnikov@lanl.gov).");
      Exceptions::amanzi_throw(msg);
    }

    Mff_cells_[c] = Mff;


    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("Matrix_MFD: unexpected failure of LAPACK in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) nokay_++;
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_PASSED) npassed_++;
  }

  // sum up the numbers across processors
  int nokay_tmp = nokay_, npassed_tmp = npassed_;
  Comm().SumAll(&nokay_tmp, &nokay_, 1);
  Comm().SumAll(&npassed_tmp, &npassed_, 1);
}


/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void
MatrixMFD::CreateMFDstiffnessMatrices(const Teuchos::Ptr<const CompositeVector>& Krel)
{
  // tag global matrices as invalid
  MarkLocalMatricesAsChanged_();

  // communicate as necessary
  if (Krel.get() && Krel->HasComponent("face")) Krel->ScatterMasterToGhosted("face");

  int dim = mesh_->getSpaceDimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;

  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  if (Aff_cells_.size() != ncells) { Aff_cells_.resize(static_cast<size_t>(ncells)); }
  if (Afc_cells_.size() != ncells) { Afc_cells_.resize(static_cast<size_t>(ncells)); }
  if (Acf_cells_.size() != ncells) { Acf_cells_.resize(static_cast<size_t>(ncells)); }
  if (Acc_cells_.size() != ncells) {
    Acc_cells_.resize(static_cast<size_t>(ncells));
    Acc_ = Teuchos::rcp(new Epetra_Vector(View, mesh_->getMap(AmanziMesh::Entity_kind::CELL,false), &Acc_cells_[0]));
  }

  for (int c = 0; c != ncells; ++c) {
    int nfaces = mesh_->getCellNumFaces(c);

    WhetStone::DenseMatrix& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (Krel == Teuchos::null || (!Krel->HasComponent("cell") && !Krel->HasComponent("face"))) {
      for (int n = 0; n != nfaces; ++n) {
        for (int m = 0; m != nfaces; ++m) { Bff(m, n) = Mff(m, n); }
      }
    } else if (Krel->HasComponent("cell") && !Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell", false);

      for (int n = 0; n != nfaces; ++n) {
        for (int m = 0; m != nfaces; ++m) { Bff(m, n) = Mff(m, n) * Krel_c[0][c]; }
      }
    } else if (!Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face", true);

      faces = mesh_->getCellFaces(c);

      for (int m = 0; m != nfaces; ++m) {
        AmanziMesh::Entity_ID f = faces[m];
        for (int n = 0; n != nfaces; ++n) { Bff(m, n) = Mff(m, n) * Krel_f[0][f]; }
      }
    } else if (Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face", true);
      const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell", false);

      faces = mesh_->getCellFaces(c);

      for (int m = 0; m != nfaces; ++m) {
        AmanziMesh::Entity_ID f = faces[m];
        for (int n = 0; n != nfaces; ++n) { Bff(m, n) = Mff(m, n) * Krel_c[0][c] * Krel_f[0][f]; }
      }
    }

    double matsum = 0.0;
    for (int n = 0; n != nfaces; ++n) {
      double rowsum = 0.0;
      double colsum = 0.0;

      for (int m = 0; m != nfaces; ++m) {
        colsum += Bff(m, n);
        rowsum += Bff(n, m);
      }

      Bcf(n) = -colsum;
      Bfc(n) = -rowsum;
      matsum += colsum;
    }

    Aff_cells_[c] = Bff;
    Afc_cells_[c] = Bfc;
    Acf_cells_[c] = Bcf;
    Acc_cells_[c] = matsum;
  }
}


/* ******************************************************************
 * Create elemental rhs vectors.
 ****************************************************************** */
void
MatrixMFD::CreateMFDrhsVectors()
{
  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  if (Ff_cells_.size() != ncells) { Ff_cells_.resize(static_cast<size_t>(ncells)); }
  if (Fc_cells_.size() != ncells) { Fc_cells_.resize(static_cast<size_t>(ncells)); }

  for (int c = 0; c != ncells; ++c) {
    int nfaces = mesh_->getCellNumFaces(c);

    Epetra_SerialDenseVector Ff(nfaces); // Entries are initilaized to 0.0.
    double Fc = 0.0;

    Ff_cells_[c] = Ff;
    Fc_cells_[c] = Fc;
  }
}


/* ******************************************************************
 * Applies boundary conditions to elemental stiffness matrices and
 * adds contributions to elemental rigth-hand-sides.
 ****************************************************************** */
void
MatrixMFD::ApplyBoundaryConditions(const std::vector<MatrixBC>& bc_markers,
                                   const std::vector<double>& bc_values,
                                   bool ADD_BC_FLUX)
{
  bc_markers_ = bc_markers;

  // tag global matrices as invalid
  MarkLocalMatricesAsChanged_();

  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);
  int nfaces = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;

  for (int c = 0; c != ncells; ++c) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];

    Epetra_SerialDenseVector& Ff = Ff_cells_[c];
    double& Fc = Fc_cells_[c];

    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      if (bc_markers[f] == MATRIX_BC_DIRICHLET) {
        for (int m = 0; m != nfaces; ++m) {
          Ff[m] -= Bff(m, n) * bc_values[f];
          Bff(n, m) = Bff(m, n) = 0.0;
        }
        Fc -= Bcf(n) * bc_values[f];

        Bcf(n) = Bfc(n) = 0.0;

        Bff(n, n) = 1.0;
        Ff[n] = bc_values[f];
      } else if ((bc_markers[f] == MATRIX_BC_FLUX) && (ADD_BC_FLUX)) {
        Ff[n] -= bc_values[f] * mesh_->getFaceArea(f);
      }
    }
  }
}


/* ******************************************************************
 * Initialize global matrices.
 *
 * This likely should only be called once.
 * If matrix is non-symmetric, we generate transpose of the matrix
 * block Afc_ to reuse cf_graph; otherwise, pointer Afc_ = Acf_.
 ****************************************************************** */
void
MatrixMFD::SymbolicAssembleGlobalMatrices()
{
  const Epetra_Map& cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL,false);
  const Epetra_Map& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,false);
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);

  int avg_entries_row = (mesh_->getSpaceDimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;

  // allocate the graphs
  Teuchos::RCP<Epetra_CrsGraph> cf_graph =
    Teuchos::rcp(new Epetra_CrsGraph(Copy, cmap, fmap_wghost, avg_entries_row, false));
  Teuchos::RCP<Epetra_FECrsGraph> ff_graph =
    Teuchos::rcp(new Epetra_FECrsGraph(Copy, fmap, 2 * avg_entries_row));

  // fill the graphs
  FillMatrixGraphs_(cf_graph.ptr(), ff_graph.ptr());

  // assemble the graphs
  int ierr = cf_graph->FillComplete(fmap, cmap);
  ASSERT(!ierr);
  ierr = ff_graph->GlobalAssemble(); // Symbolic graph is complete.
  ASSERT(!ierr);

  // allocate the matrices
  CreateMatrices_(*cf_graph, *ff_graph);
}


/* ******************************************************************
 * Fill sparsity structure graphs of global matrices (only done once)
 ****************************************************************** */
void
MatrixMFD::FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
                             const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph)
{
  const Epetra_Map& cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL,false);
  const Epetra_Map& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,false);
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);

  AmanziMesh::Entity_ID_List faces;
  int faces_LID[MFD_MAX_FACES]; // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // fill the graphs
  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);
  for (int c = 0; c != ncells; ++c) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n != nfaces; ++n) { faces_GID[n] = fmap_wghost.GID(faces[n]); }
    cf_graph->InsertMyIndices(c, nfaces, &(faces[0]));
    ff_graph->InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
  }
}


/* ******************************************************************
 * Allocate global matrices
 ****************************************************************** */
void
MatrixMFD::CreateMatrices_(const Epetra_CrsGraph& cf_graph, const Epetra_FECrsGraph& ff_graph)
{
  // create global matrices
  const Epetra_Map& cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL,false);
  Aff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Sff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Aff_->GlobalAssemble();
  Sff_->GlobalAssemble();

  // create the RHS
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::Entity_kind::CELL;
  locations[1] = AmanziMesh::Entity_kind::FACE;

  std::vector<int> num_dofs(2, 1);
  CompositeVectorSpace space;
  space.SetMesh(mesh_)->SetGhosted()->SetComponents(names, locations, num_dofs);
  rhs_ = Teuchos::rcp(new CompositeVector(space));
}


/* ******************************************************************
 * Parallel matvec product Y <-- A * X.
 ****************************************************************** */
int
MatrixMFD::Apply(const CompositeVector& X, CompositeVector& Y) const
{
  if (!Y.Ghosted()) {
    ASSERT(0);
    return 1;
  }
  if (!X.Ghosted()) {
    ASSERT(0);
    return 1;
  }

  X.ScatterMasterToGhosted();
  Y.ViewComponent("face", true)->PutScalar(0.);
  Y.ViewComponent("cell", true)->PutScalar(0.);

  const std::vector<Teuchos::SerialDenseMatrix<int, double>>& Aff = Aff_cells();
  const std::vector<Epetra_SerialDenseVector>& Afc = Afc_cells();
  const std::vector<Epetra_SerialDenseVector>& Acf = Acf_cells();

  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

  Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    Teuchos::SerialDenseVector<int, double> v(nfaces), av(nfaces);
    for (int n = 0; n < nfaces; n++) { v(n) = Xf[0][faces[n]]; }

    av.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Aff[c], v, 0.0);

    double tmp = Xc[0][c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      Yf[0][f] += av(n);
      Yc[0][c] += Acf[c](n) * v(n);
      Yf[0][f] += Afc[c](n) * tmp;
    }
    Yc[0][c] += (*Acc_)[c] * tmp;
  }
  Y.GatherGhostedToMaster("face", Add);
  return 0;
}

/* ******************************************************************
 * Parallel solve, Y <-- A^-1 X
 ****************************************************************** */
int
MatrixMFD::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  if (!assembled_schur_) {
    AssembleSchur_();
    UpdatePreconditioner_();
  }

  if (S_pc_ == Teuchos::null) {
    Errors::Message msg(
      "MatrixMFD::ApplyInverse called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }

  // Temporary cell and face vectors.
  CompositeVector T(X, true);

  // FORWARD ELIMINATION:  Tf = Xf - Afc_ inv(Acc_) Xc
  int ierr;
  Epetra_MultiVector& Tc = *T.ViewComponent("cell", false);
  ierr = Tc.ReciprocalMultiply(1.0, *Acc_, *X.ViewComponent("cell", false), 0.0);
  ASSERT(!ierr);

  ApplyAfc(T, T, 0.0);
  Epetra_MultiVector& Tf = *T.ViewComponent("face", false);
  Tf.Update(1.0, *X.ViewComponent("face", false), -1.0);

  // Solve the Schur complement system Sff_ * Yf = Tf.
  ierr = S_pc_->ApplyInverse(Tf, *Y.ViewComponent("face", false));
  ASSERT(!ierr);

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc_) (Xc - Acf_ Yf)
  ApplyAcf(Y, T, 0.0);

  Tc.Update(1.0, *X.ViewComponent("cell", false), -1.0);
  ierr |= Y.ViewComponent("cell", false)->ReciprocalMultiply(1.0, *Acc_, Tc, 0.0);

  if (ierr) {
    Errors::Message msg("MatrixMFD::ApplyInverse has failed in calculating y = A*x.");
    Exceptions::amanzi_throw(msg);
  }

  return ierr;
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = f - A * x
 ****************************************************************** */
void
MatrixMFD::ComputeResidual(const CompositeVector& solution,
                           const Teuchos::Ptr<CompositeVector>& residual) const
{
  Apply(solution, *residual);
  if (!assembled_rhs_) AssembleRHS_();
  residual->Update(1.0, *rhs_, -1.0);
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = A * x - f
 ****************************************************************** */
void
MatrixMFD::ComputeNegativeResidual(const CompositeVector& solution,
                                   const Teuchos::Ptr<CompositeVector>& residual) const
{
  Apply(solution, *residual);
  if (!assembled_rhs_) AssembleRHS_();

  residual->Update(-1.0, *rhs_, 1.0);
}


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
void
MatrixMFD::InitPreconditioner()
{}


/* ******************************************************************
 * Rebuild preconditioner.
 ****************************************************************** */
void
MatrixMFD::UpdatePreconditioner_() const
{
  if (S_pc_ == Teuchos::null) {
    Errors::Message msg(
      "MatrixMFD::ApplyInverse() called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }
  S_pc_->Destroy();

  // dump the schur complement
  // std::stringstream filename_s2;
  // filename_s2 << "schur_PC_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *Sff_);

  S_pc_->Update(Sff_);
}


/* ******************************************************************
 * WARNING: Routines requires original mass matrices (Aff_cells_), i.e.
 * before boundary conditions were imposed.
 *
 * WARNING: Since diffusive flux is not continuous, we derive it only
 * once (using flag) and in exactly the same manner as in routine
 * Flow_PK::addGravityFluxes_DarcyFlux.
 *
 ****************************************************************** */
void
MatrixMFD::DeriveFlux(const CompositeVector& solution,
                      const Teuchos::Ptr<CompositeVector>& flux) const
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;


  flux->PutScalar(0.);

  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);
  int nfaces_owned = flux->size("face", false);
  solution.ScatterMasterToGhosted("face");

  std::vector<bool> done(nfaces_owned, false);
  const Epetra_MultiVector& soln_cells = *solution.ViewComponent("cell", false);
  const Epetra_MultiVector& soln_faces = *solution.ViewComponent("face", true);
  Epetra_MultiVector& flux_v = *flux->ViewComponent("face", false);

  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->getCellFacesAndDirections(c, &faces, &dirs);
    int nfaces = faces.size();

    dp.resize(nfaces);
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      dp[n] = soln_cells[0][c] - soln_faces[0][f];
    }

    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      if (f < nfaces_owned && !done[f]) {
        double s = 0.0;
        for (int m = 0; m != nfaces; ++m) { s += Aff_cells_[c](n, m) * dp[m]; }

        flux_v[0][f] = s * dirs[n];
        done[f] = true;
      }
    }
  }

  // ensure post-condition - we got them all
  for (int f = 0; f != nfaces_owned; ++f) { ASSERT(done[f]); }
}


/* ******************************************************************
 * Derive Darcy velocity in cells.
 * WARNING: It cannot be consistent with the Darcy flux.
 * WARNING: It is assumed that flux faces have been communicated.
 ****************************************************************** */
void
MatrixMFD::DeriveCellVelocity(const CompositeVector& flux,
                              const Teuchos::Ptr<CompositeVector>& velocity) const
{
  flux.ScatterMasterToGhosted("face");
  const Epetra_MultiVector& flux_f = *flux.ViewComponent("face", true);
  Epetra_MultiVector& velocity_c = *velocity->ViewComponent("cell", false);

  int dim = mesh_->getSpaceDimension();
  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziGeometry::Point gradient(dim);
  AmanziMesh::Entity_ID_List faces;

  for (int c = 0; c != ncells_owned; ++c) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();
    std::vector<double> solution(nfaces);

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      solution[n] = flux_f[0][f];
    }

    mfd.RecoverGradient_MassMatrix(c, solution, gradient);
    for (int i = 0; i < dim; i++) velocity_c[i][c] = -gradient[i];
  }
}


/* ******************************************************************
 * Solve the bottom row of the block system for lambda, given p.
 ****************************************************************** */
void
MatrixMFD::UpdateConsistentFaceConstraints(const Teuchos::Ptr<CompositeVector>& u)
{
  if (!assembled_operator_) AssembleAff_();
  if (!assembled_rhs_) AssembleRHS_();

  // Aff solutions
  if (Aff_solver_ == Teuchos::null) {
    if (plist_.isSublist("consistent face solver")) {
      Teuchos::ParameterList Aff_plist = plist_.sublist("consistent face solver");
      Aff_op_ = Teuchos::rcp(new EpetraMatrixDefault<Epetra_FECrsMatrix>(Aff_plist));
      Aff_op_->Update(Aff_);

      if (Aff_plist.isParameter("iterative method")) {
        AmanziSolvers::LinearOperatorFactory<EpetraMatrix, Epetra_Vector, Epetra_BlockMap> op_fac;
        Aff_solver_ = op_fac.Create(Aff_plist, Aff_op_);
      } else {
        Aff_solver_ = Aff_op_;
      }
    } else {
      Errors::Message msg("MatrixMFD::UpdateConsistentFaceConstraints was called, but no "
                          "consistent face solver sublist was provided.");
      Exceptions::amanzi_throw(msg);
    }
  }

  Teuchos::Ptr<const CompositeVector> rhs = rhs_.ptr();
  const Epetra_MultiVector& rhs_f = *rhs->ViewComponent("face", false);

  Teuchos::RCP<CompositeVector> work = Teuchos::rcp(new CompositeVector(*rhs_));

  ApplyAfc(*u, *work, 0.); // Afc is kept in the transpose form.
  work->ViewComponent("face", false)->Update(1.0, rhs_f, -1.0);

  Aff_op_->Destroy();
  Aff_op_->Update(Aff_);
  int ierr = Aff_solver_->ApplyInverse(*(*work->ViewComponent("face", false))(0),
                                       *(*u->ViewComponent("face", false))(0));
}


/* ******************************************************************
 * Solve the bottom row of the block system for lambda, given p.
 ****************************************************************** */
void
MatrixMFD::UpdateConsistentFaceCorrection(const CompositeVector& u,
                                          const Teuchos::Ptr<CompositeVector>& Pu)
{
  if (!assembled_operator_) AssembleAff_();

  // Aff solutions
  if (Aff_solver_ == Teuchos::null) {
    if (plist_.isSublist("consistent face solver")) {
      Teuchos::ParameterList Aff_plist = plist_.sublist("consistent face solver");
      Aff_op_ = Teuchos::rcp(new EpetraMatrixDefault<Epetra_FECrsMatrix>(Aff_plist));
      Aff_op_->Update(Aff_);

      if (Aff_plist.isParameter("iterative method")) {
        AmanziSolvers::LinearOperatorFactory<EpetraMatrix, Epetra_Vector, Epetra_BlockMap> op_fac;
        Aff_solver_ = op_fac.Create(Aff_plist, Aff_op_);
      } else {
        Aff_solver_ = Aff_op_;
      }
    } else {
      Errors::Message msg("MatrixMFD::UpdateConsistentFaceConstraints was called, but no "
                          "consistent face solver sublist was provided.");
      Exceptions::amanzi_throw(msg);
    }
  }

  Teuchos::RCP<CompositeVector> work = Teuchos::rcp(new CompositeVector(*Pu));
  ApplyAfc(*Pu, *work, 0.); // Afc is kept in the transpose form.
  work->ViewComponent("face", false)->Update(1.0, *u.ViewComponent("face", false), -1.0);

  Aff_op_->Destroy();
  Aff_op_->Update(Aff_);
  int ierr = Aff_solver_->ApplyInverse(*(*work->ViewComponent("face", false))(0),
                                       *(*Pu->ViewComponent("face", false))(0));
  ASSERT(!ierr);
}


/* ******************************************************************
 * Solve the top row of the block system for p, given lambda.
 ****************************************************************** */
void
MatrixMFD::UpdateConsistentCellCorrection(const CompositeVector& u,
                                          const Teuchos::Ptr<CompositeVector>& Pu)
{
  Epetra_MultiVector Tc(*Pu->ViewComponent("cell", false));

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc_) (Xc - Acf_ Yf)
  ApplyAcf(*Pu, *Pu, 0.);
  Epetra_MultiVector& Pu_c = *Pu->ViewComponent("cell", false);
  Pu_c.Update(1.0, *u.ViewComponent("cell", false), -1.0);

  for (int c = 0; c != Pu_c.MyLength(); ++c) { Pu_c[0][c] /= Acc_cells_[c]; }
}


void
MatrixMFD::Add2MFDstiffnessMatrices(std::vector<double>* Acc_ptr,
                                    std::vector<Teuchos::SerialDenseMatrix<int, double>>* Aff_ptr,
                                    std::vector<Epetra_SerialDenseVector>* Acf_ptr,
                                    std::vector<Epetra_SerialDenseVector>* Afc_ptr)
{
  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  if (Acc_ptr) {
    for (int c = 0; c != ncells; ++c) { Acc_cells_[c] += (*Acc_ptr)[c]; }
  }
  if (Afc_ptr) {
    for (int c = 0; c != ncells; ++c) { Afc_cells_[c] += Afc_ptr->at(c); }
  }
  if (Acf_ptr) {
    for (int c = 0; c != ncells; ++c) { Acf_cells_[c] += Acf_ptr->at(c); }
  }
  if (Aff_ptr) {
    for (int c = 0; c != ncells; ++c) { Aff_cells_[c] += Aff_ptr->at(c); }
  }
}


/* ******************************************************************
 * Public method: Y_c = scalar * Y_c + Acf * X_f
 ****************************************************************** */
int
MatrixMFD::ApplyAcf(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  return ApplyAcf_(X, *Y.ViewComponent("cell", false), scalar);
}


int
MatrixMFD::ApplyAcf_(const CompositeVector& X, Epetra_MultiVector& Y, double scalar) const
{
  if (!X.Ghosted()) {
    ASSERT(0);
    return 1;
  }
  X.ScatterMasterToGhosted("face", true); // force scatter for now... --etc
  ApplyAcf_(*X.ViewComponent("face", true), Y, scalar);
  return 0;
}


/* ******************************************************************
 * Protected method: Y = scalar * Y + Acf * X
 ****************************************************************** */
int
MatrixMFD::ApplyAcf_(const Epetra_MultiVector& X, Epetra_MultiVector& Y, double scalar) const
{
  if (scalar == 0.0) {
    Y.PutScalar(0.0);
  } else if (scalar != 1.0) {
    Y.Scale(scalar);
  }

  const std::vector<Epetra_SerialDenseVector>& Acf = Acf_cells();

  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) { Y[0][c] += Acf[c][n] * X[0][faces[n]]; }
  }
  return 0;
}


/* ******************************************************************
 * Public method: Y_f = scalar * Y_f + Afc * X_c
 ****************************************************************** */
int
MatrixMFD::ApplyAfc(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  return ApplyAfc_(*X.ViewComponent("cell", false), Y, scalar);
}

int
MatrixMFD::ApplyAfc_(const Epetra_MultiVector& X, CompositeVector& Y, double scalar) const
{
  if (!Y.Ghosted()) {
    ASSERT(0);
    return 1;
  }

  ApplyAfc_(X, *Y.ViewComponent("face", true), scalar);
  Y.GatherGhostedToMaster("face", Add);
  return 0;
}


/* ******************************************************************
 * Protected method: Y = scalar * Y + Afc * X
 ****************************************************************** */
int
MatrixMFD::ApplyAfc_(const Epetra_MultiVector& X, Epetra_MultiVector& Y, double scalar) const
{
  if (scalar == 0.0) {
    Y.PutScalar(0.0);
  } else if (scalar != 1.0) {
    Y.Scale(scalar);

    int nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::OWNED);
    int nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::USED);
    for (int f = nfaces_owned; f < nfaces_wghost; f++) Y[0][f] = 0.0;
  }

  const std::vector<Epetra_SerialDenseVector>& Afc = Afc_cells();

  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    double tmp = X[0][c];
    for (int n = 0; n < nfaces; n++) { Y[0][faces[n]] += Afc[c][n] * tmp; }
  }
  return 0;
}


/* ******************************************************************
 * Assemble elemental rhs matrices into global RHS
 ****************************************************************** */
void
MatrixMFD::AssembleRHS_() const
{
  rhs_->ViewComponent("face", true)->PutScalar(0.0);
  Epetra_MultiVector& rhs_c = *rhs_->ViewComponent("cell", false);
  Epetra_MultiVector& rhs_f = *rhs_->ViewComponent("face", true);

  // loop over cells and fill
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
  int faces_LID[MFD_MAX_FACES];
  int faces_GID[MFD_MAX_FACES];

  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);
  for (int c = 0; c != ncells; ++c) {
    AmanziMesh::Entity_ID_List faces;
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    // assemble rhs (and simultaneously get GIDs of faces
    rhs_c[0][c] = Fc_cells_[c];
    for (int n = 0; n != nfaces; ++n) {
      AmanziMesh::Entity_ID f = faces[n];
      rhs_f[0][f] += Ff_cells_[c][n];
    }
  }

  rhs_->GatherGhostedToMaster("face");
  assembled_rhs_ = true;
}


/* ******************************************************************
 * Convert elemental mass matrices into stiffness matrices and
 * assemble them into four global matrix Aff.
 ****************************************************************** */
void
MatrixMFD::AssembleAff_() const
{
  ASSERT(Aff_.get()); // precondition: matrices have been created

  // reinitialize to zero if adding
  Aff_->PutScalar(0.0);

  AmanziMesh::Entity_ID_List faces;
  int gid[MFD_MAX_FACES];

  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  for (int c = 0; c != ncells; ++c) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n != nfaces; ++n) { gid[n] = fmap_wghost.GID(faces[n]); }
    Aff_->SumIntoGlobalValues(nfaces, gid, Aff_cells_[c].values());
  }

  // communicate
  Aff_->GlobalAssemble();

  // tag matrices as assembled
  assembled_operator_ = true;
}


/* ******************************************************************
 * Assemble Schur complement from elemental matrices.
 ****************************************************************** */
void
MatrixMFD::AssembleSchur_() const
{
  const std::vector<Teuchos::SerialDenseMatrix<int, double>>& Aff = Aff_cells();
  const std::vector<Epetra_SerialDenseVector>& Afc = Afc_cells();
  const std::vector<Epetra_SerialDenseVector>& Acf = Acf_cells();
  const std::vector<double>& Acc = Acc_cells();

  // initialize to zero
  Sff_->PutScalar(0.0);

  // loop over cells and assemble
  AmanziMesh::Entity_ID_List faces;
  const Epetra_Map& fmap_wghost = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::OWNED);

  for (int c = 0; c != ncells; ++c) {
    faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();
    Epetra_SerialDenseMatrix Tff(nfaces, nfaces); // T implies local S
    const Epetra_SerialDenseVector& Bcf = Acf[c];
    const Epetra_SerialDenseVector& Bfc = Afc[c];

    for (int n = 0; n != nfaces; ++n) {
      for (int m = 0; m != nfaces; ++m) {
        Tff(n, m) = Aff_cells_[c](n, m) - Bfc[n] * Bcf[m] / Acc[c];
      }
    }

    Epetra_IntSerialDenseVector gid(nfaces);
    for (int n = 0; n != nfaces; ++n) { // boundary conditions
      int f = faces[n];
      gid[n] = fmap_wghost.GID(f);

      if (bc_markers_[f] == MATRIX_BC_DIRICHLET) {
        for (int m = 0; m != nfaces; ++m) Tff(n, m) = Tff(m, n) = 0.0;
        Tff(n, n) = 1.0;
      }
    }

    Sff_->SumIntoGlobalValues(gid, Tff);
  }
  Sff_->GlobalAssemble();

  // tag matrices as assembled
  assembled_schur_ = true;
}

} // namespace Operators
} // namespace Amanzi
