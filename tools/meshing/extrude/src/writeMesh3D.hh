/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef MESH_WRITER_HH_
#define MESH_WRITER_HH_

#include "Mesh3D.hh"


namespace Amanzi {
namespace AmanziGeometry {

void writeMesh3D_exodus(const Mesh3D& m, const std::string& filename);

}
} // namespace Amanzi


#endif
