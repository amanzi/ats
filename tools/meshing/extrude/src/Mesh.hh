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


struct Mesh3D {
  Mesh3D(const Mesh2D * const m_, int n_layers);

  void extrude(double dz, const std::vector<int>& cell_set) {
    std::vector<double> dzs(m->coords.size(), dz);
    extrude(dzs, cell_set);
  }

  void extrude(const std::vector<double>& dzs, int my_cell_set) {
    std::vector<int> cell_set(m->ncells, my_cell_set);
    extrude(dzs, cell_set);
  }

  void extrude(double dz, int my_cell_set) {
    std::vector<double> dzs(m->coords.size(), dz);
    std::vector<int> cell_set(m->ncells, my_cell_set);
    extrude(dzs, cell_set);
  }

  int node_structure(int n_2d, int layer) {
    return n_2d + layer*m->coords.size();
}

  int cell_structure(int c_2d, int layer) {
    return c_2d + layer*m->cell2node.size();
  }

  void extrude(const std::vector<double>& dz,
               const std::vector<int>& cell_set);
  void finish_sets(); 

  const Mesh2D * const m;

  std::vector<Point> coords;
  std::vector<std::vector<int> > cell2face;
  std::vector<std::vector<int> > face2node;

  std::vector<int> block_ids;
  
  std::vector<int> side_face_counts;

  std::vector<std::pair<std::vector<int>,
                        std::vector<int> > > face_sets;
  std::vector<int> face_sets_id;


  int current_layer;
  int total_layers;
  Point datum;
};

}
}

#endif
