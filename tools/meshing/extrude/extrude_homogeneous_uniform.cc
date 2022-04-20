#include "Mesh3D.hh"
#include "writeMesh3D.hh"
#include "readMesh2D.hh"


int main() {
  using namespace Amanzi::AmanziGeometry;

  std::string mesh_in = "Mesh.txt";
  std::string mesh_out = "Mesh3D_Homogeneous2mSoil.exo";
  
  std::vector<double> ref_soil_mlay_dz = {0.02, 0.03, 0.05, 0.15, 0.25, 0.5, 1.0};
  std::vector<double> ref_bedrock_mlay_dz = {10., 10.};

  int hmg_soil_type = 1000;
  int hmg_bedrock_type = 100;
  
  int nsoil_lay = ref_soil_mlay_dz.size();
  int nbedrock_lay = ref_bedrock_mlay_dz.size();

  std::cout << "Extruding: " << mesh_in << " and writing to: " << mesh_out << std::endl;

  std::vector<int> soil_type;
  std::vector<int> bedrock_type;
  std::vector<double> depths;

  //  auto m = readMesh2D_text(mesh_in, soil_type, bedrock_type, depths);
  auto m = readMesh2D_text(mesh_in, soil_type, bedrock_type, depths, 630497, 4.04083e+06);

  int nsnodes = m.coords.size();
  
  Mesh3D m3(&m, nsoil_lay + nbedrock_lay);
  for (int ilay=0; ilay!=nsoil_lay; ++ilay) {
    m3.extrude(ref_soil_mlay_dz[ilay], hmg_soil_type);
  }
  for (int ilay=0; ilay!=nbedrock_lay; ++ilay) {
    m3.extrude(ref_bedrock_mlay_dz[ilay], hmg_bedrock_type);
  }
  
  m3.finish();

  std::cout << "NNodes on the surf = " << m.coords.size() << std::endl;
  std::cout << "Ncells on the surf = " << m.cell2node.size() << std::endl;
  std::cout << "NNodes on 3D = " << m3.coords.size() << std::endl;
  std::cout << "Ncells on 3D = " << m3.cell2face.size() << std::endl;

  writeMesh3D_exodus(m3, mesh_out);
  return 0;
}
