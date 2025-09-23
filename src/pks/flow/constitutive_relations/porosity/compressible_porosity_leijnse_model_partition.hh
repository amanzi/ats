/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  A collection of porosity models along with a Mesh Partition.
*/

#ifndef AMANZI_FLOW_RELATIONS_COMP_PORO_LEIJNSE_PARTITION_
#define AMANZI_FLOW_RELATIONS_COMP_PORO_LEIJNSE_PARTITION_

#include "compressible_porosity_leijnse_model.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

typedef std::vector<Teuchos::RCP<CompressiblePorosityLeijnseModel>>
  CompressiblePorosityLeijnseModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, CompressiblePorosityLeijnseModelList>
  CompressiblePorosityLeijnseModelPartition;

// Non-member factory
Teuchos::RCP<CompressiblePorosityLeijnseModelPartition>
createCompressiblePorosityLeijnseModelPartition(Teuchos::ParameterList& plist);

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
