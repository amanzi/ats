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
                mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED));
  int d = mesh.space_dimension() - 1;
  for (int col = 0; col != mesh.num_columns(); ++col) {
    const auto& col_faces = mesh.faces_of_column(col);
    const auto& col_cells = mesh.cells_of_column(col);
    double top_z = mesh.face_centroid(col_faces[0])[d];
    for (int i = 0; i != col_cells.size(); ++i) {
      double cell_z =
        (mesh.face_centroid(col_faces[i])[d] + mesh.face_centroid(col_faces[i + 1])[d]) / 2;
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
                mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED));
  int d = mesh.space_dimension() - 1;
  for (int col = 0; col != mesh.num_columns(); ++col) {
    const auto& col_faces = mesh.faces_of_column(col);
    const auto& col_cells = mesh.cells_of_column(col);
    double top_z = mesh.face_centroid(col_faces[0])[d];
    for (int i = 0; i != col_cells.size(); ++i) {
      double cell_z = mesh.cell_centroid(col_cells[i])[d];
      depth[0][col_cells[i]] = top_z - cell_z;
    }
  }

  double mv;
  depth.MinValue(&mv);
  AMANZI_ASSERT(mv > 0.);
}

} // namespace Flow
} // namespace Amanzi
