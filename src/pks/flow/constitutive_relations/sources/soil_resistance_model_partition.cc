/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*
  Mesh Partition for soil type determined soil resistance.
  Currently this is only used for Sakagucki-Zeng model, but
  probably useful for other soil type determined soil resistance
  models that might be added in future.
*/


#include "dbc.hh"
#include "soil_resistance_model_partition.hh"


namespace Amanzi {
namespace Flow {

// Non-member factory
Teuchos::RCP<SoilResistanceModelPartition>
createSoilResistanceModelPartition(Teuchos::ParameterList& plist)
{
  SoilResistanceSakaguckiZengModelList rs_list;
  std::vector<std::string> region_list;

  for (Teuchos::ParameterList::ConstIterator lcv = plist.begin(); lcv != plist.end(); ++lcv) {
    std::string name = lcv->first;
    if (plist.isSublist(name)) {
      Teuchos::ParameterList sublist = plist.sublist(name);
      region_list.push_back(sublist.get<std::string>("region"));
      rs_list.push_back(Teuchos::rcp(new SoilResistanceSakaguckiZengModel(sublist)));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  Teuchos::RCP<Functions::MeshPartition> part =
    Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::CELL, region_list));

  return Teuchos::rcp(new SoilResistanceModelPartition(part, rs_list));
}

} // namespace Flow
} // namespace Amanzi
