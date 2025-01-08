/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  A collection of porosity models along with a Mesh Partition.
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
