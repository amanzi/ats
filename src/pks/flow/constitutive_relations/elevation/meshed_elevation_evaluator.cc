/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  An elevation evaluator getting values from the volumetric mesh.

*/
#define _USE_MATH_DEFINES
#include <cmath>

#include "MeshAlgorithms.hh"
#include "Mesh.hh"
#include "Point.hh"
#include "meshed_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

namespace Impl {

void
slope_aspect(const AmanziGeometry::Point& normal, double& slope, double& aspect)
{
  // -- S = || n - (n dot z) z || / | n dot z |
  slope = std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2)) / std::abs(normal[2]);

  // and aspect
  if (normal[0] > 0.0) {
    // right half
    if (normal[1] > 0.0) {
      // upper right quadrant
      aspect = std::atan(normal[0] / normal[1]);
    } else if (normal[1] < 0.0) {
      // lower right quadrant
      aspect = M_PI - std::atan(normal[0] / -normal[1]);
    } else {
      // due east
      aspect = M_PI_2;
    }
  } else if (normal[0] < 0.0) {
    // left half
    if (normal[1] > 0.0) {
      // upper left quadrant
      aspect = 2 * M_PI - std::atan(-normal[0] / normal[1]);
    } else if (normal[1] < 0.0) {
      // lower left quadrant
      aspect = M_PI + std::atan(normal[0] / -normal[1]);
    } else {
      // due west
      aspect = 3 * M_PI_2;
    }
  } else {
    // north or south
    if (normal[1] > 0.0) {
      aspect = 0.0;
    } else if (normal[1] < 0.0) {
      aspect = M_PI;
    } else {
      // normal is (0,0,1)
      aspect = 0.;
    }
  }
}

} // namespace Impl


MeshedElevationEvaluator::MeshedElevationEvaluator(Teuchos::ParameterList& plist)
  : ElevationEvaluator(plist){};


Teuchos::RCP<Evaluator>
MeshedElevationEvaluator::Clone() const
{
  return Teuchos::rcp(new MeshedElevationEvaluator(*this));
}


void
MeshedElevationEvaluator::EvaluateElevationAndSlope_(const State& S,
                                                     const std::vector<CompositeVector*>& results)
{
  CompositeVector* elev = results[0];
  CompositeVector* slope = results[1];
  CompositeVector* aspect = results[2];

  Epetra_MultiVector& elev_c = *elev->ViewComponent("cell", false);
  Epetra_MultiVector& slope_c = *slope->ViewComponent("cell", false);
  Epetra_MultiVector& aspect_c = *aspect->ViewComponent("cell", false);

  // Get the elevation and slope values from the domain mesh.
  Key domain = Keys::getDomain(my_keys_.front().first);
  const auto& surface_mesh = S.GetMesh(domain);
  const auto& parent_mesh = surface_mesh->getParentMesh();
  AMANZI_ASSERT(parent_mesh != Teuchos::null);
  AMANZI_ASSERT(parent_mesh->getSpaceDimension() == 3);

  if (parent_mesh->getManifoldDimension() == 3) {
    // surface mesh is extraction
    int ncells = elev_c.MyLength();
    for (int c = 0; c != ncells; ++c) {
      // Set the elevation on cells by getting the corresponding face and its
      // centroid.
      auto domain_face = surface_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
      elev_c[0][c] = parent_mesh->getFaceCentroid(domain_face)[2];

      // Set the slope and aspect using the upward normal of the corresponding
      // face.
      AmanziGeometry::Point normal = parent_mesh->getFaceNormal(domain_face);
      Impl::slope_aspect(normal, slope_c[0][c], aspect_c[0][c]);
    }

    // Set the elevation on faces by getting the corresponding nodes and
    // averaging.
    if (elev->HasComponent("face")) {
      Epetra_MultiVector& elev_f = *elev->ViewComponent("face", false);
      AmanziMesh::Entity_ID node0, node1;
      AmanziGeometry::Point coord0(3);
      AmanziGeometry::Point coord1(3);

      int nfaces = elev_f.MyLength();
      for (int f = 0; f != nfaces; ++f) {
        auto surface_nodes = surface_mesh->getFaceNodes(f);
        node0 = surface_mesh->getEntityParent(AmanziMesh::Entity_kind::NODE, surface_nodes[0]);
        node1 = surface_mesh->getEntityParent(AmanziMesh::Entity_kind::NODE, surface_nodes[1]);
        coord0 = parent_mesh->getNodeCoordinate(node0);
        coord1 = parent_mesh->getNodeCoordinate(node1);
        elev_f[0][f] = (coord0[2] + coord1[2]) / 2.0;
      }
    }

    if (elev->HasComponent("boundary_face")) {
      const Epetra_Map& vandalay_map =
        surface_mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false);
      const Epetra_Map& face_map = surface_mesh->getMap(AmanziMesh::Entity_kind::FACE, false);
      Epetra_MultiVector& elev_bf = *elev->ViewComponent("boundary_face", false);

      AmanziMesh::Entity_ID node0, node1;
      AmanziGeometry::Point coord0(3);
      AmanziGeometry::Point coord1(3);

      for (int bf = 0; bf != elev_bf.MyLength(); ++bf) {
        AmanziMesh::Entity_ID f = AmanziMesh::getBoundaryFaceFace(*surface_mesh, bf);
        auto surface_nodes = surface_mesh->getFaceNodes(f);
        node0 = surface_mesh->getEntityParent(AmanziMesh::Entity_kind::NODE, surface_nodes[0]);
        node1 = surface_mesh->getEntityParent(AmanziMesh::Entity_kind::NODE, surface_nodes[1]);
        coord0 = parent_mesh->getNodeCoordinate(node0);
        coord1 = parent_mesh->getNodeCoordinate(node1);
        elev_bf[0][bf] = (coord0[2] + coord1[2]) / 2.0;
      }
    }

  } else {
    // surface mesh is a flattening
    int ncells = elev_c.MyLength();
    for (int c = 0; c != ncells; ++c) {
      // Set the elevation on cells by getting the corresponding face and its
      // centroid.
      auto domain_cell = surface_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
      elev_c[0][c] = parent_mesh->getCellCentroid(domain_cell)[2];

      // Figure out the upward normal using edges.
      // -- faces of the cell
      AmanziGeometry::Point normal(0., 0., 0.);
      auto faces = parent_mesh->getCellFaces(domain_cell);

      // -- Get the normals of all faces of the surface cell.
      int count = faces.size();
      std::vector<AmanziGeometry::Point> normals(count);
      for (int lcv = 0; lcv != count; ++lcv) {
        normals[lcv] = parent_mesh->getFaceNormal(faces[lcv], domain_cell);
      }

      // -- Average the cross product of successive faces to get a cell normal.
      for (int lcv = 0; lcv != (count - 1); ++lcv) {
        normal[0] += normals[lcv][1] * normals[lcv + 1][2] - normals[lcv][2] * normals[lcv + 1][1];
        normal[1] += normals[lcv][2] * normals[lcv + 1][0] - normals[lcv][0] * normals[lcv + 1][2];
        normal[2] += normals[lcv][0] * normals[lcv + 1][1] - normals[lcv][1] * normals[lcv + 1][0];
      }
      normal[0] += normals[count - 1][1] * normals[0][2] - normals[count - 1][2] * normals[0][1];
      normal[1] += normals[count - 1][2] * normals[0][0] - normals[count - 1][0] * normals[0][2];
      normal[2] += normals[count - 1][0] * normals[0][1] - normals[count - 1][1] * normals[0][0];
      normal /= count;

      Impl::slope_aspect(normal, slope_c[0][c], aspect_c[0][c]);
    }

    // Set the elevation on faces by getting the corresponding face and its
    // centroid.
    if (elev->HasComponent("face")) {
      Epetra_MultiVector& elev_f = *elev->ViewComponent("face", false);
      int nfaces = elev_f.MyLength();
      for (int f = 0; f != nfaces; ++f) {
        // Note that a surface face is a surface mesh's face.
        AmanziMesh::Entity_ID domain_face =
          surface_mesh->getEntityParent(AmanziMesh::Entity_kind::FACE, f);
        AmanziGeometry::Point x = parent_mesh->getFaceCentroid(domain_face);
        elev_f[0][f] = x[2];
      }
    }
  }
}

} // namespace Flow
} // namespace Amanzi
