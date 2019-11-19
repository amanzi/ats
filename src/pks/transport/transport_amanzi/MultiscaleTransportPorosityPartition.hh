/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  A collection of multiscale porosity models along with a mesh partition.
*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_PARTITION_HH_
#define MULTISCALE_TRANSPORT_POROSITY_PARTITION_HH_

#include "Mesh.hh"
#include "MeshPartition.hh"
#include "MultiscaleTransportPorosity.hh"

namespace Amanzi {
namespace Transport {

typedef std::vector<Teuchos::RCP<MultiscaleTransportPorosity> > MsTList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, MsTList> MultiscaleTransportPorosityPartition;

// Non-member factory takes the global trnsport list
Teuchos::RCP<MultiscaleTransportPorosityPartition> CreateMultiscaleTransportPorosityPartition(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    Teuchos::RCP<Teuchos::ParameterList> tplist);

// Non-member function quering partition
int NumberMatrixNodes(Teuchos::RCP<MultiscaleTransportPorosityPartition>& msp);

}  // namespace Transport
}  // namespace Amanzi

#endif
