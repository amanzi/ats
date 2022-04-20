#include "Mesh3D.hh"
#include "writeMesh3D.hh"
#include "readMesh2D.hh"


int main() {
  using namespace Amanzi::AmanziGeometry;

  std::string mesh_in = "Mesh.txt";
  std::string mesh_out = "Mesh3D_OneLayer.exo";
  
  std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out << std::endl;

  std::vector<int> soil_type;
  std::vector<int> bedrock_type;
  std::vector<double> depths;

  auto m = readMesh2D_text(mesh_in, soil_type, bedrock_type, depths);
  int nsnodes = m.coords.size();
  
  Mesh3D m3(&m, 1);
  m3.extrude(0.02, 100);
  m3.finish();

  std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
  std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
  std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
  std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;

  writeMesh3D_exodus(m3, mesh_out);
  return 0;
}
