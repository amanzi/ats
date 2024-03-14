/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A collection of WRMs along with a Mesh Partition.
/*!

A WRM partition is a list of (region, WRM) pairs, where the regions partition
the mesh.

.. _wrm-partition-typedinline-spec:
.. admonition:: wrm-partition-typedinline-spec

   * `"region`" ``[string]`` Region on which the WRM is valid.
   * `"WRM type`" ``[string]`` Name of the WRM type.
   * `"_WRM_type_ parameters`" ``[_WRM_type_-spec]`` See below for the required
     parameter spec for each type.

*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PARTITION_
#define AMANZI_FLOW_RELATIONS_WRM_PARTITION_

#include "wrm.hh"
#include "wrm_permafrost_model.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<WRM>> WRMList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, WRMList> WRMPartition;

typedef std::vector<Teuchos::RCP<WRMPermafrostModel>> WRMPermafrostModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, WRMPermafrostModelList>
  WRMPermafrostModelPartition;

// Non-member factory
Teuchos::RCP<WRMPartition>
createWRMPartition(Teuchos::ParameterList& plist);

Teuchos::RCP<WRMPermafrostModelPartition>
createWRMPermafrostModelPartition(Teuchos::ParameterList& plist, Teuchos::RCP<WRMPartition>& wrms);

} // namespace Flow
} // namespace Amanzi

#endif
