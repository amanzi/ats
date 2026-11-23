/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

/*
  Provides a collection of models along with a Mesh Partition: 
  per region with mesh partition of that region.
*/

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Flow {

template<class Model>
using RegionModelList = std::vector<Teuchos::RCP<Model>>;

template<class Model>
using RegionModelPartition =
  std::pair<Teuchos::RCP<Functions::MeshPartition>, RegionModelList<Model>>;


// Non-member factory
template<class Model>
Teuchos::RCP<RegionModelPartition<Model>>
createRegionModelPartition(Teuchos::ParameterList& plist)
{
  RegionModelList<Model> model_list;
  std::vector<std::string> region_list;

  for (Teuchos::ParameterList::ConstIterator lcv = plist.begin(); lcv != plist.end(); ++lcv) {
    std::string name = lcv->first;
    if (plist.isSublist(name)) {
      Teuchos::ParameterList sublist = plist.sublist(name);
      region_list.push_back(sublist.get<std::string>("region"));
      model_list.push_back(Teuchos::rcp(new Model(sublist)));
    } else {
      AMANZI_ASSERT(0);
    }
  }

  Teuchos::RCP<Functions::MeshPartition> part =
    Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::CELL, region_list));
  return Teuchos::rcp(new RegionModelPartition<Model>(part, model_list));
}

} // namespace Flow
} // namespace Amanzi
