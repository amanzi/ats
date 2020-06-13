/*
  Evaluates depth of various mesh entities.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_DEPTH_MODEL_HH_
#define AMANZI_FLOW_DEPTH_MODEL_HH_


#include "Mesh.hh"

namespace Amanzi {
namespace Flow {

void
DepthModel(const AmanziMesh::Mesh& mesh, Epetra_MultiVector& depth);

void
DepthModel_Cell(int c, const AmanziMesh::Mesh& mesh,
                Epetra_MultiVector& depth);

} //namespace
} //namespace

#endif
