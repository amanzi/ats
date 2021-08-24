#include <set>
#include "Mesh2D.hh"

namespace Amanzi {
namespace AmanziGeometry {

Mesh2D::Mesh2D(std::vector<Point>& coords_,
               std::vector<std::vector<int> >& cell2node_,
               std::vector<std::vector<int> >& cell_sets_) :
    coords(coords_),
    cell2node(cell2node_),
    cell_sets(cell_sets_),
    datum(2)
{
  double area_eps = 1.e-6;
  
  for (auto& set : cell_sets) {
    AMANZI_ASSERT(set.size() == cell2node.size());
  }

  nnodes = coords.size();
  ncells = cell2node.size();
  
  for (auto& c : cell2node) {
    Point v1(2), v2(2);
    v1.set(coords[c[1]][0] - coords[c[0]][0], coords[c[1]][1] - coords[c[0]][1]);
    v2.set(coords[c[2]][0] - coords[c[0]][0], coords[c[2]][1] - coords[c[0]][1]);
    Point cross = v1^v2;
    if (is_greater(0.0, cross[0])) {
      std::reverse(c.begin(), c.end());
    }
    if (is_greater(area_eps, std::abs(cross[0]))) {
      std::cout << "Zero area triangle:" << std::endl
                << " " << coords[c[0]] << std::endl
                << " " << coords[c[1]] << std::endl
                << " " << coords[c[2]] << std::endl;
      AMANZI_ASSERT(false);
    }
  }

  for (int c=0; c!=cell2node.size(); ++c) {
    std::vector<int> c2f(3);
    
    std::vector<int> f1 = {cell2node[c][0], cell2node[c][1]};
    c2f[0] = face_constructor(f1, 0, c);
    std::vector<int> f2 = {cell2node[c][1], cell2node[c][2]};
    c2f[1] = face_constructor(f2, 1, c);
    std::vector<int> f3 = {cell2node[c][2], cell2node[c][0]};
    c2f[2] = face_constructor(f3, 2, c);      
    cell2face.emplace_back(c2f);
  }

  // set sizes
  nfaces = face2node.size();
  AMANZI_ASSERT(ncells == cell2face.size());

  // set boundaries
  std::vector<int> boundary_c, boundary_f;
  for (int lcv=0; lcv!=side_face_counts.size(); ++lcv) {
    if (side_face_counts[lcv] == 1) {
      boundary_c.push_back(face_cell_when_created[lcv]);
      boundary_f.push_back(face_in_cell_when_created[lcv]);
    }
  }
  boundary_faces = std::make_pair(boundary_c,boundary_f);

  // normalize the nodes
  int64_t x=0., y=0.;
  for (auto& p : coords) {
    x += p[0];
    y += p[1];
  }
  datum[0] = x / coords.size();
  datum[1] = y / coords.size();
  for (auto& p : coords) {
    p[0] -= datum[0];
    p[1] -= datum[1];
  }
}

int
Mesh2D::face_constructor(const std::vector<int>& nodes,
                         int index_in_cell,
                         int cell)
{
  auto h = nodes[0] > nodes[1] ? hash(nodes[1], nodes[0]) :
      hash(nodes[0],nodes[1]);

  auto match = faces_sorted.find(h);
  if (match != faces_sorted.end()) {
    int f = match->second;
    side_face_counts[f]++;
    return f;
  }

  // not already created
  int f = face2node.size();
  faces_sorted[h] = f;
  face2node.emplace_back(nodes);
  face_in_cell_when_created.push_back(index_in_cell);
  face_cell_when_created.push_back(cell);
  side_face_counts.push_back(1);
  return f;
}



}
}
