#include "Mesh.hh"
#include "writeMesh3D.hh"
#include "readMesh2D.hh"


int main() {
  using namespace Amanzi::AmanziGeometry;

  std::string mesh_in = "Mesh.txt";
  std::string mesh_out = "Mesh3D.exo";

  std::vector<int> soil_type;
  std::vector<int> bedrock_type;
  std::vector<int> veg_type;
  std::vector<double> depths;

  Mesh2D m = readMesh2D_text(mesh_in, soil_type, veg_type, bedrock_type, depths);

  Mesh3D m3(m, 3);
  m3.extrude(depths, soil_type);
  m3.extrude(10, 1000);
  m3.extrude(10, 1000);
  
  writeMesh3D_exodus(m3, mesh_out);
  return 0;
}
