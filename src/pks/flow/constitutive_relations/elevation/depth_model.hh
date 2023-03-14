/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates depth of various mesh entities.

*/

#ifndef AMANZI_FLOW_DEPTH_MODEL_HH_
#define AMANZI_FLOW_DEPTH_MODEL_HH_

#include "Epetra_MultiVector.h"

#include "Mesh.hh"

namespace Amanzi {
namespace Flow {

void
DepthModel(const AmanziMesh::Mesh& mesh, Epetra_MultiVector& depth);

void
DepthModel_Cell(int c, const AmanziMesh::Mesh& mesh, Epetra_MultiVector& depth);

} // namespace Flow
} // namespace Amanzi

#endif
