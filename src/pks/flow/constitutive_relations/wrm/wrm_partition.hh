/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PARTITION_
#define AMANZI_FLOW_RELATIONS_WRM_PARTITION_

#include "wrm.hh"
#include "wrm_permafrost_model.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace ATS_Physics {      
namespace Flow {

typedef std::vector<Teuchos::RCP<WRM>> WRMList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, WRMList> WRMPartition;

typedef std::vector<Teuchos::RCP<WRMPermafrostModel>> WRMPermafrostModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, WRMPermafrostModelList>
  WRMPermafrostModelPartition;

// Non-member factory
Teuchos::RCP<WRMPartition> createWRMPartition(Teuchos::ParameterList& plist);

Teuchos::RCP<WRMPermafrostModelPartition> createWRMPermafrostModelPartition(
  Teuchos::ParameterList& plist,
  Teuchos::RCP<WRMPartition>& wrms);

} // namespace Flow
}
} // namespace Amanzi

#endif
