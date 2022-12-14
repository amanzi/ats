/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  A collection of porosity models along with a Mesh Partition.
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_COMP_PORO_LEIJNSE_PARTITION_
#define AMANZI_FLOW_RELATIONS_COMP_PORO_LEIJNSE_PARTITION_

#include "compressible_porosity_leijnse_model.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<CompressiblePorosityLeijnseModel>>
  CompressiblePorosityLeijnseModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, CompressiblePorosityLeijnseModelList>
  CompressiblePorosityLeijnseModelPartition;

// Non-member factory
Teuchos::RCP<CompressiblePorosityLeijnseModelPartition>
createCompressiblePorosityLeijnseModelPartition(Teuchos::ParameterList& plist);

} // namespace Flow
} // namespace Amanzi

#endif
