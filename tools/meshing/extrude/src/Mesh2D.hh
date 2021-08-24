#ifndef MESH_2D_HH_
#define MESH_2D_HH_

#include <vector>
#include <utility>
#include <algorithm>
#include <map>

#include "Point.hh"
#include <stdint.h>

namespace Amanzi {
namespace AmanziGeometry {

struct Mesh2D {
  Mesh2D(std::vector<Point>& coords_,
         std::vector<std::vector<int> >& cell2node_,
         std::vector<std::vector<int> >& cell_sets_);
  
  int face_constructor(const std::vector<int>& nodes,
                       int face_in_cell,
                       int cell);

  int64_t hash(int i, int j) {
    return ((int64_t) (nnodes+1))*((int64_t) i) + (int64_t) j;
  }
  
  std::vector<Point> coords;
  std::vector<std::vector<int> > cell2node;
  std::vector<std::vector<int> > cell2face;
  std::vector<std::vector<int> > face2node;

  std::vector<std::vector<int> > cell_sets;
  std::pair<std::vector<int>,
            std::vector<int> > boundary_faces;
  
  std::map<int64_t, int> faces_sorted;
  std::vector<int> side_face_counts;
  std::vector<int> face_in_cell_when_created;
  std::vector<int> face_cell_when_created;

  int nnodes, ncells, nfaces;

  Point datum;
};


}
}

#endif
