/// Used for testing 0 dz extrusion only!
// Extrudes 1 layer, with dz = 0 in one corner.
//
#include <cassert>
#include "Mesh3D.hh"
#include "writeMesh3D.hh"
#include "readMesh2D.hh"


int main() {
  using namespace Amanzi::AmanziGeometry;

  std::string mesh_in = "Mesh.txt";

  std::vector<int> soil_type;
  std::vector<int> bedrock_type;
  std::vector<double> depths;

  //  auto m = readMesh2D_text(mesh_in, soil_type, bedrock_type, depths);
  auto m = readMesh2D_text(mesh_in, soil_type, bedrock_type, depths, 630497, 4.04083e+06);
  int nsnodes = m.coords.size();
  assert(nsnodes == 4);
  assert(m.cell2node.size() == 2);

  {  
    std::string mesh_out = "Mesh3D_Pinchout.exo";
    std::string mesh_out_ns = "Mesh3D_Pinchout_NonSquashed.exo";
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out << std::endl;

    std::vector<double> dzs = { 1., 0., 0.5, 1. }; 
  
    // make an extruded with pinchouts, squashed
    Mesh3D m3(&m, 1);
    m3.extrude(1.0, 1001, true);
    m3.extrude(dzs, 1001, true);
    m3.finish();

    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;
    assert(m3.coords.size() == 11);
    assert(m3.cell2face.size() == 4);
    assert(m3.face2node.size() == 16);
    
    writeMesh3D_exodus(m3, mesh_out);

    // also a non-squashed version
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out_ns << std::endl;
    Mesh3D m3_ns(&m, 1);
    m3_ns.extrude(1.0, 1001, false);
    m3_ns.extrude(dzs, 1001, false);
    m3_ns.finish();
    
    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3_ns.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3_ns.cell2face.size() << std::endl;
    assert(m3_ns.coords.size() == 12);
    assert(m3_ns.cell2face.size() == 4);
    assert(m3_ns.face2node.size() == 16);
    
    writeMesh3D_exodus(m3_ns, mesh_out_ns);
  }



  {  
    std::string mesh_out = "Mesh3D_Pinchout2.exo";
    std::string mesh_out_ns = "Mesh3D_Pinchout2_NonSquashed.exo";
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out << std::endl;

    std::vector<double> dzs = { 0., 0., 0.5, 1. }; 
  
    // make an extruded with pinchouts, squashed
    Mesh3D m3(&m, 1);
    m3.extrude(1.0, 1001, true);
    m3.extrude(dzs, 1001, true);
    m3.finish();

    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;
    assert(m3.coords.size() == 10);
    assert(m3.cell2face.size() == 4);
    assert(m3.face2node.size() == 15);
    
    writeMesh3D_exodus(m3, mesh_out);

    // also a non-squashed version
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out_ns << std::endl;
    Mesh3D m3_ns(&m, 1);
    m3_ns.extrude(1.0, 1001, false);
    m3_ns.extrude(dzs, 1001, false);
    m3_ns.finish();
    
    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3_ns.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3_ns.cell2face.size() << std::endl;
    assert(m3_ns.coords.size() == 12);
    assert(m3_ns.cell2face.size() == 4);
    assert(m3_ns.face2node.size() == 16);
    
    writeMesh3D_exodus(m3_ns, mesh_out_ns);
  }


  {  
    std::string mesh_out = "Mesh3D_Pinchout3.exo";
    std::string mesh_out_ns = "Mesh3D_Pinchout3_NonSquashed.exo";
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out << std::endl;

    std::vector<double> dzs = { 0., 0., 0., 1. }; 
  
    // make an extruded with pinchouts, squashed
    Mesh3D m3(&m, 1);
    m3.extrude(1.0, 1001, true);
    m3.extrude(dzs, 1001, true);
    m3.finish();

    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;
    assert(m3.coords.size() == 9);
    assert(m3.cell2face.size() == 3);
    assert(m3.face2node.size() == 12);
    
    writeMesh3D_exodus(m3, mesh_out);

    // also a non-squashed version
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out_ns << std::endl;
    Mesh3D m3_ns(&m, 1);
    m3_ns.extrude(1.0, 1001, false);
    m3_ns.extrude(dzs, 1001, false);
    m3_ns.finish();
    
    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3_ns.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3_ns.cell2face.size() << std::endl;
    assert(m3_ns.coords.size() == 12);
    assert(m3_ns.cell2face.size() == 4);
    assert(m3_ns.face2node.size() == 16);
    
    writeMesh3D_exodus(m3_ns, mesh_out_ns);
  }

  {  
    std::string mesh_out = "Mesh3D_Pinchout4.exo";
    std::string mesh_out_ns = "Mesh3D_Pinchout4_NonSquashed.exo";
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out << std::endl;

    std::vector<double> dzs = { 0., 0., 0., 0. }; 
  
    // make an extruded with pinchouts, squashed
    Mesh3D m3(&m, 1);
    m3.extrude(1.0, 1001, true);
    m3.extrude(dzs, 1001, true);
    m3.finish();

    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;
    assert(m3.coords.size() == 8);
    assert(m3.cell2face.size() == 2);
    assert(m3.face2node.size() == 9);
    
    writeMesh3D_exodus(m3, mesh_out);

    // also a non-squashed version
    std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out_ns << std::endl;
    Mesh3D m3_ns(&m, 1);
    m3_ns.extrude(1.0, 1001, false);
    m3_ns.extrude(dzs, 1001, false);
    m3_ns.finish();
    
    std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
    std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
    std::cout << "NNodes on 3D = " << m3_ns.coords.size() << std::endl;
    std::cout << "Ncells on 3D = " << m3_ns.cell2face.size() << std::endl;
    assert(m3_ns.coords.size() == 12);
    assert(m3_ns.cell2face.size() == 4);
    assert(m3_ns.face2node.size() == 16);
    
    writeMesh3D_exodus(m3_ns, mesh_out_ns);
  }
  
  return 0;
}
