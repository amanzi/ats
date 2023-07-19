/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov) (ATS version)
*/

/*
  License:

  Deformation optimization matrix
*/

#include <set>
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "CompositeVectorSpace.hh"
#include "MatrixVolumetricDeformation.hh"
#include "PreconditionerFactory.hh"

#define MESH_TYPE 1 // 0 = HEXES, 1 = TRIANGULAR PRISMS

namespace Amanzi {
namespace Operators {


MatrixVolumetricDeformation::MatrixVolumetricDeformation(
  Teuchos::ParameterList& plist,
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : plist_(plist), mesh_(mesh)
{
  InitializeFromOptions_();
  PreAssemble_();
};


MatrixVolumetricDeformation::MatrixVolumetricDeformation(const MatrixVolumetricDeformation& other)
  : plist_(other.plist_), mesh_(other.mesh_)
{
  InitializeFromOptions_();
  PreAssemble_();
};


// Apply matrix, b <-- Ax
int
MatrixVolumetricDeformation::Apply(const CompositeVector& x, CompositeVector& b) const
{
  return operator_->Apply(*x.ViewComponent("node", false), *b.ViewComponent("node", false));
}


// Apply the inverse, x <-- A^-1 b
int
MatrixVolumetricDeformation::ApplyInverse(const CompositeVector& b, CompositeVector& x) const
{
  int ierr = prec_->ApplyInverse(*b.ViewComponent("node", false), *x.ViewComponent("node", false));
  AMANZI_ASSERT(!ierr);
  return ierr;
}


void
MatrixVolumetricDeformation::InitializeFromOptions_()
{
  // parameters for optimization
  smoothing_ = plist_.get<double>("smoothing coefficient");
  diagonal_shift_ = plist_.get<double>("diagonal shift", 1.e-6);

  // preconditioner
  AmanziPreconditioners::PreconditionerFactory fac;
  prec_ = fac.Create(plist_);
};


// This is a Normal equation, so we need to apply N^T to the rhs
void
MatrixVolumetricDeformation::ApplyRHS(
  const CompositeVector& x_cell,
  const Teuchos::Ptr<CompositeVector>& x_node,
  const Teuchos::Ptr<const AmanziMesh::Entity_ID_List>& fixed_nodes) const
{
  // Equations of the form: dVdz x = b, solved via:
  //   dVdz^T * dVdz * x = dVdz^T * b, where operator_ = dVdz^T * dVdz
  // -- form dVdz^T * b
  dVdz_->SetUseTranspose(true);
  dVdz_->Apply(*x_cell.ViewComponent("cell", false), *x_node->ViewComponent("node", false));

  // Must also apply the boundary condition, dz_bottom = 0
  // -- Fix the bottom nodes, they may not move
  Epetra_MultiVector& rhs_n = *x_node->ViewComponent("node", false);
  unsigned int nnodes_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  for (AmanziMesh::Entity_ID_List::const_iterator n = fixed_nodes->begin(); n != fixed_nodes->end();
       ++n) {
    if (*n < nnodes_owned) rhs_n[0][*n] = 0;
  }
}


void
MatrixVolumetricDeformation::PreAssemble_()
{
  int ierr = 0;

  // HARD CODED CRUFT!
#if MESH_TYPE
  int nnz = 6;
  int indices[6];
  double values[6];
#else
  int nnz = 8;
  int indices[8];
  double values[8];
#endif


  double eps = 1.e-4;

  // Domain and Range: matrix inverse solves the problem, given dV, what is
  // dnode_z that results in that dV.  Therefore, A * dnode_z = dV
  const Epetra_Map& cell_map = mesh_->getMap(AmanziMesh::Entity_kind::CELL,false);
  const Epetra_Map& node_map = mesh_->getMap(AmanziMesh::Entity_kind::NODE,false);
  const Epetra_Map& node_map_wghost = mesh_->getMap(AmanziMesh::Entity_kind::NODE,true);

  range_ = Teuchos::rcp(new CompositeVectorSpace());
  range_->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  domain_ = Teuchos::rcp(new CompositeVectorSpace());
  domain_->SetMesh(mesh_)->SetComponent("node", AmanziMesh::Entity_kind::NODE, 1);

  unsigned int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  unsigned int nnodes = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);


  // == Create dVdz ==
  // -- create the matrix
  dVdz_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cell_map, nnz, true));

  // -- Assemble
  for (unsigned int c = 0; c != ncells; ++c) {
    for (int n = 0; n != nnz; ++n) values[n] = 0.;

    AmanziMesh::Entity_ID_List nodes;
    nodes = mesh_->getCellNodes(c);
    AMANZI_ASSERT(nodes.size() == nnz);

    // determine the upward/downward faces
    AmanziMesh::Entity_ID_List faces;
    faces = mesh_->getCellFaces(c);

    int my_up_n = -1;
    int my_down_n = -1;
    for (int n = 0; n != faces.size(); ++n) {
      if (mesh_->getFaceNormal(faces[n], c)[2] > eps) {
        AMANZI_ASSERT(my_up_n < 0);
        my_up_n = n;
      } else if (mesh_->getFaceNormal(faces[n], c)[2] < -eps) {
        AMANZI_ASSERT(my_down_n < 0);
        my_down_n = n;
      }
    }
    AMANZI_ASSERT(my_up_n >= 0);
    AMANZI_ASSERT(my_down_n >= 0);

    // Determine the perpendicular area (remember, getFaceNormal() is weighted
    // by face area already).
    double perp_area_up = std::abs(mesh_->getFaceNormal(faces[my_up_n])[2]);
    double perp_area_down = std::abs(mesh_->getFaceNormal(faces[my_down_n])[2]);
    AMANZI_ASSERT(std::abs(perp_area_up - perp_area_down) < 1.e-6);

    // loop over nodes in the top face, setting the Jacobian
    AmanziMesh::Entity_ID_List nodes_up;
    nodes_up = mesh_->getFaceNodes(faces[my_up_n]);
    for (int n = 0; n != nodes_up.size(); ++n) {
      values[n] = perp_area_up / (nnz / 2.);
      indices[n] = node_map_wghost.GID(nodes_up[n]);
    }

    // loop over nodes in the bottom face, setting the Jacobian
    AmanziMesh::Entity_ID_List nodes_down;
    nodes_down = mesh_->getFaceNodes(faces[my_down_n]);
    int nnz_half = nnz / 2;
    for (int n = 0; n != nodes_down.size(); ++n) {
      values[n + nnz_half] = -perp_area_up / (nnz / 2.);
      indices[n + nnz_half] = node_map_wghost.GID(nodes_down[n]);
    }

    // ensure all are nonzero
    for (int n = 0; n != nnz; ++n) AMANZI_ASSERT(std::abs(values[n]) > 0.);

    // Assemble
    int ierr = dVdz_->InsertGlobalValues(cell_map.GID(c), nnz, values, indices);
    AMANZI_ASSERT(!ierr);
  }

  ierr = dVdz_->FillComplete(node_map, cell_map);
  AMANZI_ASSERT(!ierr);

  // dump dvdz
  //  EpetraExt::RowMatrixToMatlabFile("dvdz.txt", *dVdz_);

  // == Form the normal equations, dVdz^T * dVdz ==

  // calculate nnz
  int* nnz_op = new int[nnodes];
  max_nnode_neighbors_ = 0;
  for (unsigned int no = 0; no != nnodes; ++no) {
    AmanziMesh::Entity_ID_List cells;
    cells = mesh_->getNodeCells(no, AmanziMesh::Parallel_kind::ALL);
    // THIS IS VERY SPECIFIC TO GEOMETRY!
#if MESH_TYPE
    int nnode_neighbors = (cells.size() / 2 + 1) * 3;
#else
    int nnode_neighbors = 9 * 3;
#endif
    max_nnode_neighbors_ = std::max(nnode_neighbors, max_nnode_neighbors_);
    nnz_op[no] = nnode_neighbors;
  }

  operatorPre_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, node_map, nnz_op));
  delete[] nnz_op;

  // form explicitly dVdz^T * dVdz
  EpetraExt::MatrixMatrix::Multiply(*dVdz_, true, *dVdz_, false, *operatorPre_);

  // dump Operator intermediate step
  //  EpetraExt::RowMatrixToMatlabFile("op_normal_eq.txt", *operatorPre_);

  // == Add in optimization components ==
  // -- Add a small shift to the diagonal
  Epetra_Vector diag(node_map);
  ierr = operatorPre_->ExtractDiagonalCopy(diag);
  AMANZI_ASSERT(!ierr);
  for (int n = 0; n != nnodes; ++n) diag[n] += diagonal_shift_;
  ierr = operatorPre_->ReplaceDiagonalValues(diag);
  AMANZI_ASSERT(!ierr);

  // dump Operator intermediate step
  //  EpetraExt::RowMatrixToMatlabFile("op_diag.txt", *operatorPre_);

  // -- Add in a diffusive term
  int* node_indices = new int[max_nnode_neighbors_];
  double* node_values = new double[max_nnode_neighbors_];

  for (int n = 0; n != nnodes; ++n) {
    // determine the set of node neighbors
    std::set<int> neighbors;
    AmanziMesh::Entity_ID_List faces;
    mesh_->node_get_faces(n, AmanziMesh::Parallel_kind::ALL, &faces);
    for (AmanziMesh::Entity_ID_List::const_iterator f = faces.begin(); f != faces.end(); ++f) {
      AmanziMesh::Entity_ID_List neighbor_nodes;
      neighbor_nodes = mesh_->getFaceNodes(*f);
      neighbors.insert(neighbor_nodes.begin(), neighbor_nodes.end());
    }

    // remove my node
    neighbors.erase(n);

    // apply the stencil
    int nneighbors = neighbors.size() + 1;
    AmanziGeometry::Point center(3);
    center = mesh_->getNodeCoordinate(n);
    node_indices[0] = node_map.GID(n);
    node_values[0] = 0.;

    int lcv = 1;
    for (std::set<int>::const_iterator hn = neighbors.begin(); hn != neighbors.end(); ++hn) {
      node_indices[lcv] = node_map_wghost.GID(*hn);

      AmanziGeometry::Point coord(3);
      coord = mesh_->getNodeCoordinate(*hn);
      double smoothing = smoothing_ / AmanziGeometry::norm(center - coord);
      node_values[0] += smoothing;
      node_values[lcv] = -smoothing;
      lcv++;
    }

    // add into the row
    ierr =
      operatorPre_->SumIntoGlobalValues(node_map.GID(n), nneighbors, node_values, node_indices);
    AMANZI_ASSERT(!ierr);
  }

  ierr = operatorPre_->FillComplete();
  AMANZI_ASSERT(!ierr);

  // clean up
  delete[] node_indices;
  delete[] node_values;
}


void
MatrixVolumetricDeformation::Assemble(
  const Teuchos::Ptr<const AmanziMesh::Entity_ID_List>& fixed_nodes)
{
  int ierr = 0;

  // HARD CODED CRUFT!
#if MESH_TYPE
  int nnz = 6;
  int indices[6];
  double values[6];
#else
  int nnz = 8;
  int indices[8];
  double values[8];
#endif

  double eps = 1.e-4;
  int* node_indices = new int[max_nnode_neighbors_];
  double* node_values = new double[max_nnode_neighbors_];

  // Domain and Range: matrix inverse solves the problem, given dV, what is
  // dnode_z that results in that dV.  Therefore, A * dnode_z = dV
  const Epetra_Map& node_map = mesh_->getMap(AmanziMesh::Entity_kind::NODE,false);
  const Epetra_Map& node_map_wghost = mesh_->getMap(AmanziMesh::Entity_kind::NODE,true);
  unsigned int nnodes = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);

  // reset to non-fixed operator.
  if (operator_ == Teuchos::null) {
    operator_ = Teuchos::rcp(new Epetra_CrsMatrix(*operatorPre_));
  } else {
    *operator_ = *operatorPre_;
  }

  Epetra_Vector diag(node_map);
  ierr = operatorPre_->ExtractDiagonalCopy(diag);
  AMANZI_ASSERT(!ierr);
  double min(0.);
  diag.MinValue(&min);
  std::cout << "MIN VAL OP_PRE = " << min << std::endl;
  std::cout << "SIZE nnodes = " << nnodes << std::endl;
  std::cout << "SIZE map = " << node_map.NumMyElements() << std::endl;
  std::cout << "SIZE op pre Row = " << operatorPre_->RowMap().NumMyElements() << std::endl;
  std::cout << "SIZE op pre Domain = " << operatorPre_->DomainMap().NumMyElements() << std::endl;
  std::cout << "SIZE op pre Range = " << operatorPre_->RangeMap().NumMyElements() << std::endl;
  std::cout << "SIZE op Row = " << operator_->RowMap().NumMyElements() << std::endl;
  std::cout << "SIZE op Domain = " << operator_->DomainMap().NumMyElements() << std::endl;
  std::cout << "SIZE op Range = " << operator_->RangeMap().NumMyElements() << std::endl;

  ierr = operator_->ExtractDiagonalCopy(diag);
  AMANZI_ASSERT(!ierr);
  diag.MinValue(&min);
  std::cout << "MIN VAL OP (PRE) = " << min << std::endl;

  // -- Fix the bottom nodes, they may not move
  for (AmanziMesh::Entity_ID_List::const_iterator n = fixed_nodes->begin(); n != fixed_nodes->end();
       ++n) {
    // extract the row
    if (*n < nnodes) {
      int n_gid = node_map.GID(*n);
      int nneighbors;
      ierr = operator_->ExtractGlobalRowCopy(
        n_gid, max_nnode_neighbors_, nneighbors, node_values, node_indices);
      AMANZI_ASSERT(!ierr);

      // zero the row, 1 on diagonal
      for (int i = 0; i != nneighbors; ++i) node_values[i] = node_indices[i] == n_gid ? 1. : 0.;

      double rowsum = 0.;
      for (int i = 0; i != nneighbors; ++i) rowsum += node_indices[i];
      AMANZI_ASSERT(rowsum > 0.);

      // replace the row
      ierr = operator_->ReplaceGlobalValues(n_gid, nneighbors, node_values, node_indices);
      AMANZI_ASSERT(!ierr);
    }
  }

  ierr = operator_->FillComplete();
  AMANZI_ASSERT(!ierr);

  ierr = operator_->ExtractDiagonalCopy(diag);
  AMANZI_ASSERT(!ierr);
  diag.MinValue(&min);
  std::cout << "MIN VAL OP (POST) = " << min << std::endl;

  // dump Operator intermediate step
  //  EpetraExt::RowMatrixToMatlabFile("op.txt", *operator_);

  // clean up
  delete[] node_indices;
  delete[] node_values;

  // Set the operator in the precon
  prec_->Destroy();
  prec_->Update(operator_);
}


void
MatrixVolumetricDeformation::set_inverse_parameters(){};

} // namespace Operators
} // namespace Amanzi
