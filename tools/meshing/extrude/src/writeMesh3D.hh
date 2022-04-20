#ifndef MESH_WRITER_HH_
#define MESH_WRITER_HH_

#include "Mesh3D.hh"


namespace Amanzi {
namespace AmanziGeometry {

void writeMesh3D_exodus(const Mesh3D& m, const std::string& filename);

}
}


#endif
