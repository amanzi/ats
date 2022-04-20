#ifndef MESH_3D_HH_
#define MESH_3D_HH_

#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <stdint.h>

#include "Point.hh"
#include "Mesh2D.hh"


namespace Amanzi {
namespace AmanziGeometry {


struct Mesh3D {
  Mesh3D(const Mesh2D * const m_, int n_layers);

  void extrude(double dz, const std::vector<int>& cell_set, bool squash_zero_edges=true) {
    std::vector<double> dzs(m->coords.size(), dz);
    extrude(dzs, cell_set, squash_zero_edges);
  }

  void extrude(const std::vector<double>& dzs, int my_cell_set, bool squash_zero_edges=true) {
    std::vector<int> cell_set(m->ncells, my_cell_set);
    extrude(dzs, cell_set, squash_zero_edges);
  }

  void extrude(double dz, int my_cell_set, bool squash_zero_edges=true) {
    std::vector<double> dzs(m->coords.size(), dz);
    std::vector<int> cell_set(m->ncells, my_cell_set);
    extrude(dzs, cell_set, squash_zero_edges);
  }

  void extrude(const std::vector<double>& dz,
               const std::vector<int>& cell_set,
               bool squash_zero_edges=true);
  void finish(); 

  const Mesh2D * const m;

  // basic geometric/topology info
  std::vector<Point> coords;
  std::vector<std::vector<int> > cell2face;
  std::vector<std::vector<int> > face2node;

  // labels
  std::vector<int> block_ids;
  std::vector<std::pair<std::vector<int>,
                        std::vector<int> > > side_sets;
  std::vector<int> side_sets_id;

  // internally used to extrude
  std::vector<int> up_faces;
  std::vector<int> dn_faces;
  std::vector<int> up_nodes;
  std::vector<int> dn_nodes;
  std::vector<int> cells_in_col;

  // other meta-data
  int current_layer;
  int total_layers;
  Point datum;
};

}
}

#endif
