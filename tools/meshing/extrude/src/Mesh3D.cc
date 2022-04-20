#include <set>
#include <algorithm>
#include <numeric>
#include <fstream>

#include "Mesh2D.hh"
#include "Mesh3D.hh"

#define logically_structured

namespace Amanzi {
namespace AmanziGeometry {

Mesh3D::Mesh3D(const Mesh2D * const m_, int n_layers) :
    m(m_),
    current_layer(0),
    total_layers(n_layers),
    datum(m_->datum),
    cells_in_col(m->ncells, 0)    
{
  // reserve/allocate space
  Point d(3);
  int n_nodes = m->nnodes*(n_layers+1);
  coords.reserve(n_nodes);
  
  int n_cells = n_layers * m->ncells;
  cell2face.reserve(n_cells);

  int n_faces = n_layers * m->nfaces
      + (n_layers+1)*m->ncells;
  face2node.reserve(n_faces);

  // copy the top surface coords
  coords.insert(coords.end(), m->coords.begin(), m->coords.end());

  // create the top layer of faces
  face2node.insert(face2node.end(), m->cell2node.begin(),
                   m->cell2node.end());
  up_faces.resize(face2node.size());
  std::iota(up_faces.begin(), up_faces.end(), 0);
  up_nodes.resize(coords.size());
  std::iota(up_nodes.begin(), up_nodes.end(), 0);
  dn_nodes = up_nodes;
  dn_faces = up_faces;

  // create the "bottom" sideset
  side_sets.emplace_back(std::piecewise_construct,
                         std::forward_as_tuple(m->ncells, -1),
                         std::forward_as_tuple(m->ncells, 1));
  side_sets_id.push_back(1);
  
  // create the "surface" sideset
  side_sets.emplace_back(std::piecewise_construct,
                         std::forward_as_tuple(m->ncells, -1),
                         std::forward_as_tuple(m->ncells, 0));
  side_sets_id.push_back(2);

  // create an empty "sides" sideset
  side_sets.emplace_back(std::make_pair(std::vector<int>(), std::vector<int>()));
  side_sets_id.push_back(3);
}


void
Mesh3D::extrude(const std::vector<double>& dz,
                const std::vector<int>& block_ids_,
                bool squash_zero_edges) {
  AMANZI_ASSERT(dz.size() == m->coords.size());
  AMANZI_ASSERT(block_ids_.size() == m->cell2node.size());

  auto this_layer_sides = std::vector<int>(m->face2node.size(), -1);
  auto node_differs = [this](int n) {return this->dn_nodes[n] != this->up_nodes[n];};
  auto node_same_horiz = [this](int n) {return is_equal(coords[this->dn_nodes[n]][0],
                                                        coords[this->up_nodes[n]][0])
                                            && is_equal(coords[this->dn_nodes[n]][1],
                                                        coords[this->up_nodes[n]][1]);};
  // shift the up-node coordinates by dz
  for (int n=0; n!=dz.size(); ++n) {
    if (!squash_zero_edges || dz[n] > 0.) {
      Point nc(coords[up_nodes[n]]);
      nc[2] -= dz[n];
      coords.emplace_back(std::move(nc));
      dn_nodes[n] = coords.size()-1;
    }
  }

  // add cells, faces
  for (int c=0; c!=m->ncells; ++c) {
    if (std::any_of(m->cell2node[c].begin(), m->cell2node[c].end(), node_differs)) {
      cells_in_col[c]++;
      
      // add the bottom face
      std::vector<int> dn_face;
      for (auto n : m->cell2node[c]) dn_face.emplace_back(dn_nodes[n]);
      int my_dn_f = face2node.size();
      face2node.emplace_back(std::move(dn_face));
      dn_faces[c] = my_dn_f;

      // push back a cell containing the up, dn faces
      auto cell_faces = std::vector<int>{up_faces[c], dn_faces[c]};
      int my_c = cell2face.size();

      // if this is the top cell, put it into the surface side set
      if (side_sets[1].first[c] < 0) side_sets[1].first[c] = my_c;
      // put this cell into the bottom side set -- will be overwritten if any lower
      side_sets[0].first[c] = my_c;

      // add faces for the sides as needed
      for (auto sf : m->cell2face[c]) {
        int my_f = this_layer_sides[sf];

        AMANZI_ASSERT(node_same_horiz(m->face2node[sf][0]));
        AMANZI_ASSERT(node_same_horiz(m->face2node[sf][1]));
        
        if (std::any_of(m->face2node[sf].begin(), m->face2node[sf].end(),
                        node_differs)) {
          if (my_f < 0) {
            // may need to create the face
            my_f = face2node.size();
            auto side_nodes = std::vector<int>{ up_nodes[m->face2node[sf][1]], 
                                                up_nodes[m->face2node[sf][0]] };
            if (node_differs(m->face2node[sf][0])) side_nodes.push_back(dn_nodes[m->face2node[sf][0]]);
            if (node_differs(m->face2node[sf][1])) side_nodes.push_back(dn_nodes[m->face2node[sf][1]]);
            face2node.emplace_back(std::move(side_nodes));
            cell_faces.push_back(my_f);
            this_layer_sides[sf] = my_f;

            // check if this is a boundary side, and add it to the side_set if so
            if (m->side_face_counts[sf] == 1) {
              side_sets[2].first.push_back(my_c);
              side_sets[2].second.push_back(cell_faces.size() - 1);
            }
            
          } else {
            // no need to create the face, but the cell still needs to know it
            cell_faces.push_back(my_f);
          }

        }
      }

      // finally add the cell, including a block id
      cell2face.emplace_back(std::move(cell_faces));
      block_ids.push_back(block_ids_[c]);
    }
  }

  // increment the layer metadata
  current_layer++;
  up_nodes = dn_nodes;
  up_faces = dn_faces;

  AMANZI_ASSERT(block_ids.size() == cell2face.size());
  std::cout << "POST-Extruding: currently " << cell2face.size() << " cells and " << face2node.size() << " faces." << std::endl;

}

void
Mesh3D::finish() {
  // flip the bottom faces for proper outward orientation
  for (auto f : dn_faces)
    std::reverse(face2node[f].begin(), face2node[f].end());

  // move the 2d cell sets to face sets on the surface
  std::set<int> set_ids;
  for (auto& part : m->cell_sets) {
    set_ids.insert(part.begin(), part.end());
  }

  for (int sid : set_ids) {
    std::vector<int> set_cells;
    for (auto& part : m->cell_sets) {
      for (int c=0; c!=part.size(); ++c) {
        if (part[c] == sid) {
          set_cells.push_back(side_sets[1].first[c]);
        }
      }
    }
    std::vector<int> set_faces(set_cells.size(), 0);
    side_sets.emplace_back(std::make_pair(std::move(set_cells),
            std::move(set_faces)));
    side_sets_id.push_back(sid);
  }

  // check side sets
  std::vector<int> side_face_counts(face2node.size(), 0);
  for (auto& c : cell2face)
    for (auto& f : c)
      side_face_counts[f]++;
  
  for (int lcv_s=0; lcv_s!=side_sets.size(); ++lcv_s) {
    auto& fs = side_sets[lcv_s]; 
    for (int i=0; i!=side_sets[lcv_s].first.size(); ++i) {
      int c = fs.first[i];
      int fi = fs.second[i];
      int f = cell2face[c][fi];
      if (side_face_counts[f] != 1) {
        std::cout << "Face Set " << side_sets_id[lcv_s] << ": face = " << f << " (" << c << "," << fi << ") has been counted " << side_face_counts[f] << " times (should be 1)!" << std::endl;
      }
    }
  }

  std::ofstream fid;
  fid.open("col_counts.txt");
  for (auto c : cells_in_col) fid << c << std::endl;
  fid.close();
}


}
}
