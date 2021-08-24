#include <set>
#include <vector>
#include <algorithm>
#include "exodusII.h"

#include "dbc.hh"

#include "writeMesh3D.hh"


namespace Amanzi {
namespace AmanziGeometry {

void
writeMesh3D_exodus(const Mesh3D& m, const std::string& filename) {
  // create the exodus file
  int CPU_word_size = sizeof(float);
  int IO_Word_size = 8;
  int fid = ex_create(filename.c_str(), EX_NOCLOBBER, &CPU_word_size, &IO_Word_size);
  if (fid < 0) {
    std::cerr << "Cowardly not clobbering: \"" << filename << "\" already exists." << std::endl;
    return;
  }
  
  // make the blocks by set
  std::set<int> set_ids(m.block_ids.begin(), m.block_ids.end());
  std::vector<std::vector<int> > blocks;
  std::vector<int> blocks_ncells;
  std::vector<int> blocks_id;
  std::vector<std::vector<int> > block_face_counts;

  int new_id = 0;
  std::vector<int> cell_map(m.cell2face.size(), -1);
  for (auto sid : set_ids) {
    std::vector<int> block;
    std::vector<int> face_counts;
    // count things
    int ncells = 0;

    for (int i=0; i!=m.cell2face.size(); ++i) {
      if (m.block_ids[i] == sid) {
        ncells++;
        block.insert(block.end(), m.cell2face[i].begin(), m.cell2face[i].end());
        face_counts.push_back(m.cell2face[i].size());
        cell_map[i] = new_id;
        new_id++;
      }
    }

    blocks.emplace_back(block);
    blocks_ncells.push_back(ncells);
    blocks_id.push_back(sid);
    block_face_counts.emplace_back(face_counts);
  }


  // create the params
  ex_init_params params;
  sprintf(params.title, "my_mesh");
  params.num_dim = 3;
  params.num_nodes = m.coords.size();
  params.num_edge = 0;
  params.num_edge_blk = 0;
  params.num_face = m.face2node.size();
  params.num_face_blk = 1;
  params.num_elem = m.cell2face.size();
  params.num_elem_blk = set_ids.size();
  params.num_node_maps = 0;
  params.num_edge_maps = 0;
  params.num_face_maps = 0;
  params.num_elem_maps = 0;
  params.num_side_sets = m.side_sets.size();
  params.num_elem_sets = 0;
  params.num_node_sets = 0;
  params.num_face_sets = 0;
  params.num_edge_sets = 0;

  int ierr = ex_put_init_ext(fid, &params);
  AMANZI_ASSERT(!ierr);
  
  // make the coordinate arrays, set the coordinates
  // NOTE: exodus seems to only deal with floats!
  std::vector<std::vector<float> > coords(3);
  for (int i=0; i!=3; ++i) {
    coords[i].resize(m.coords.size());
  }

  for (int n=0; n!=coords[0].size(); ++n) {
    coords[0][n] = m.coords[n][0];
    coords[1][n] = m.coords[n][1];
    coords[2][n] = m.coords[n][2];
  }

  char* coord_names[3];
  char a[10]="xcoord";
  char b[10]="ycoord";
  char c[10]="zcoord";
  coord_names[0]=a;
  coord_names[1]=b;
  coord_names[2]=c;

  ierr |= ex_put_coord_names(fid, coord_names);
  AMANZI_ASSERT(!ierr);
  ierr |= ex_put_coord(fid, &coords[0][0], &coords[1][0], &coords[2][0]);
  AMANZI_ASSERT(!ierr);

  
  // put in the face block
  std::vector<int> facenodes;
  std::vector<int> facenodes_counts;
  for (auto& nodes : m.face2node) {
    facenodes_counts.push_back(nodes.size());
    facenodes.insert(facenodes.end(), nodes.begin(), nodes.end());
  }
  ierr |= ex_put_block(fid, EX_FACE_BLOCK, 1, "NSIDED",
                       m.face2node.size(), facenodes.size(), 0,0,0);
  AMANZI_ASSERT(!ierr);

  ierr |= ex_put_entity_count_per_polyhedra(fid, EX_FACE_BLOCK, 1,
          &facenodes_counts[0]);
  AMANZI_ASSERT(!ierr);

  for (auto& e : facenodes) e++;
  ierr |= ex_put_conn(fid, EX_FACE_BLOCK, 1, &facenodes[0], NULL, NULL);
  AMANZI_ASSERT(!ierr);
  

  // put in the element blocks
  for (int lcvb=0; lcvb!=blocks.size(); ++lcvb) {
    ierr |= ex_put_block(fid, EX_ELEM_BLOCK, blocks_id[lcvb], "NFACED",
                         blocks_ncells[lcvb], 0, 0, blocks[lcvb].size(),0);
    AMANZI_ASSERT(!ierr);

    ierr |= ex_put_entity_count_per_polyhedra(fid, EX_ELEM_BLOCK, blocks_id[lcvb],
            &block_face_counts[lcvb][0]);
    AMANZI_ASSERT(!ierr);

    for (auto&e : blocks[lcvb]) e++;
    ierr |= ex_put_conn(fid, EX_ELEM_BLOCK, blocks_id[lcvb], NULL, NULL, &blocks[lcvb][0]);
    AMANZI_ASSERT(!ierr);
  }

  
  // add the side sets, mapping elems to the new ids
  for (int lcvs=0; lcvs!=m.side_sets.size(); ++lcvs) {
    auto& s = m.side_sets[lcvs];
    std::vector<int> elems_copy(s.first.size(), -1);
    auto faces_copy(s.second);
    for (int i=0; i!=elems_copy.size(); ++i) {
      elems_copy[i] = cell_map[s.first[i]] + 1;
    }
    for (auto& e : faces_copy) e++;
    ierr |= ex_put_set_param(fid, EX_SIDE_SET, m.side_sets_id[lcvs], elems_copy.size(), 0);
    AMANZI_ASSERT(!ierr);
    ierr |= ex_put_set(fid, EX_SIDE_SET, m.side_sets_id[lcvs], &elems_copy[0], &faces_copy[0]);
    AMANZI_ASSERT(!ierr);
  }

  ierr |= ex_close(fid);
  AMANZI_ASSERT(!ierr);


  // debugging/nice output
  std::cout << "Wrote 3D Mesh:" << std::endl
            << "  ncells = " << m.cell2face.size() << std::endl
            << "  nfaces = " << m.face2node.size() << std::endl
            << "  nnodes = " << m.coords.size() << std::endl
            << std::endl
            << "  side sets = " << std::endl;
  for (int i=0; i!=m.side_sets.size(); ++i)
    std::cout << "    " << m.side_sets_id[i] << " ("
              << m.side_sets[i].first.size() << " faces)" << std::endl;
  std::cout << std::endl
            << "  block ids = " << std::endl;
  for (int i=0; i!=blocks_id.size(); ++i)
    std::cout << "    " << blocks_id[i] << " ("
              << blocks_ncells[i] << " cells)" << std::endl;
  std::cout << std::endl;
  
}

}
}
