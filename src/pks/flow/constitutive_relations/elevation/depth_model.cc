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
DepthModel(const AmanziMesh::Mesh& mesh, Epetra_MultiVector& depth)
{
  depth.PutScalar(-1);
  AMANZI_ASSERT(depth.MyLength() ==
                mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED));
  for (int c = 0; c != depth.MyLength(); ++c) {
    if (depth[0][c] <= 0.) { DepthModel_Cell(c, mesh, depth); }
  }
}


void
DepthModel_Cell(int c, const AmanziMesh::Mesh& mesh, Epetra_MultiVector& depth)
{
  int z_dim = mesh.space_dimension() - 1;
  int c_above = mesh.cell_get_cell_above(c);
  if (c_above == -1) {
    // top cell, find the face above
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh.cell_get_faces_and_dirs(c, &faces, &dirs);
    int f_above = -1;
    for (auto f : faces) {
      AmanziGeometry::Point face_normal = mesh.face_normal(f, false, c);
      face_normal /= AmanziGeometry::norm(face_normal);
      if (face_normal[z_dim] > 1.e-10) {
        f_above = f;
        break;
      }
    }

    // get the depth
    depth[0][c] = mesh.face_centroid(f_above)[z_dim] - mesh.cell_centroid(c)[z_dim];
    return;
  }

  if (depth[0][c_above] <= 0) { DepthModel_Cell(c_above, mesh, depth); }
  AMANZI_ASSERT(depth[0][c_above] > 0.);
  depth[0][c] =
    depth[0][c_above] + mesh.cell_centroid(c_above)[z_dim] - mesh.cell_centroid(c)[z_dim];
  return;
}


} // namespace Flow
} // namespace Amanzi
