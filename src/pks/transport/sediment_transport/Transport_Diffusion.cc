/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "MFD3D_Diffusion.hh"
//#include "nlfv.hh"
#include "Tensor.hh"
//#include "PreconditionerFactory.hh"

#include "sediment_transport_pk.hh"
#include "SedimentTransportDefs.hh"

namespace Amanzi {
namespace SedimentTransport {

/* *******************************************************************
* Calculate dispersive tensor from given Darcy fluxes. The flux is
* assumed to be scaled by face area.
******************************************************************* */
void
SedimentTransport_PK::CalculateDiffusionTensor_(const Epetra_MultiVector& km,
                                                const Epetra_MultiVector& ws,
                                                const Epetra_MultiVector& mol_density)
{
  D_.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);

  //AmanziGeometry::Point velocity(dim);
  AmanziMesh::Entity_ID_List faces;
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  for (int c = 0; c < ncells_owned; ++c) {
    double mol_den = mol_density[0][c];
    double ponded_depth = ws[0][c];
    D_[c].PutScalar(km[0][c] * mol_den * ponded_depth);
  }
}


} // namespace SedimentTransport
} // namespace Amanzi
