/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  A collection of comp poro models along with a Mesh Partition.
*/

#include "dbc.hh"
#include "compressible_porosity_model_partition.hh"


namespace Amanzi {
namespace Flow {

// Non-member factory
Teuchos::RCP<CompressiblePorosityModelPartition>
createCompressiblePorosityModelPartition(Teuchos::ParameterList& plist)
{
  CompressiblePorosityModelList mlist;
  std::vector<std::string> region_list;

  for (Teuchos::ParameterList::ConstIterator lcv = plist.begin(); lcv != plist.end(); ++lcv) {
    std::string name = lcv->first;
    if (plist.isSublist(name)) {
      Teuchos::ParameterList sublist = plist.sublist(name);
      region_list.push_back(sublist.get<std::string>("region"));
      mlist.push_back(Teuchos::rcp(new CompressiblePorosityModel(sublist)));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  Teuchos::RCP<Functions::MeshPartition> part =
    Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::Entity_kind::CELL, region_list));

  return Teuchos::rcp(new CompressiblePorosityModelPartition(part, mlist));
}

} // namespace Flow
} // namespace Amanzi
