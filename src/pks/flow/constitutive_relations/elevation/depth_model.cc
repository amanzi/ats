/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates depth of various mesh entities.

*/

#include "depth_model.hh"

namespace Amanzi {
namespace Flow {

void
computeDepth_MeanFaceCentroid(const AmanziMesh::Mesh& mesh, Epetra_MultiVector& depth)
{
  depth.PutScalar(-1);
  AMANZI_ASSERT(depth.MyLength() ==
                mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
  int d = mesh.getSpaceDimension() - 1;
  for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
    const auto& col_faces = mesh.columns.getFaces(col);
    const auto& col_cells = mesh.columns.getCells(col);
    double top_z = mesh.getFaceCentroid(col_faces[0])[d];
    for (int i = 0; i != col_cells.size(); ++i) {
      double cell_z =
        (mesh.getFaceCentroid(col_faces[i])[d] + mesh.getFaceCentroid(col_faces[i + 1])[d]) / 2;
      depth[0][col_cells[i]] = top_z - cell_z;
    }
  }

  double mv;
  depth.MinValue(&mv);
  AMANZI_ASSERT(mv > 0.);
}


void
computeDepth_CellCentroid(const AmanziMesh::Mesh& mesh, Epetra_MultiVector& depth)
{
  depth.PutScalar(-1);
  AMANZI_ASSERT(depth.MyLength() ==
                mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
  int d = mesh.getSpaceDimension() - 1;
  for (int col = 0; col != mesh.columns.num_columns_owned; ++col) {
    const auto& col_faces = mesh.columns.getFaces(col);
    const auto& col_cells = mesh.columns.getCells(col);
    double top_z = mesh.getFaceCentroid(col_faces[0])[d];
    for (int i = 0; i != col_cells.size(); ++i) {
      double cell_z = mesh.getCellCentroid(col_cells[i])[d];
      depth[0][col_cells[i]] = top_z - cell_z;
    }
  }

  double mv;
  depth.MinValue(&mv);
  AMANZI_ASSERT(mv > 0.);
}

} // namespace Flow
} // namespace Amanzi
