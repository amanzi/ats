/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  A collection of porosity models along with a Mesh Partition.
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_COMP_PORO_PARTITION_
#define AMANZI_FLOW_RELATIONS_COMP_PORO_PARTITION_

#include "compressible_porosity_model.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<CompressiblePorosityModel>> CompressiblePorosityModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, CompressiblePorosityModelList>
  CompressiblePorosityModelPartition;

// Non-member factory
Teuchos::RCP<CompressiblePorosityModelPartition>
createCompressiblePorosityModelPartition(Teuchos::ParameterList& plist);

} // namespace Flow
} // namespace Amanzi

#endif
