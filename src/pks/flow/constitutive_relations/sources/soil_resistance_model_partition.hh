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
*/

/*!

The partition is a list of (region, WRM parameters) pairs, where the
regions partition the mesh. Note that the (region, WRM parameters) pairs
are not required arguments in the soil resistance evaluator, but
defined through model parameters under state to ensure consistency
of (region, WRM parameters) pairs for WRM and/or relative permeability
uses through the whole input file.

*/


#ifndef AMANZI_FLOW_RELATIONS_SOIL_RESISTANCE_PARTITION_
#define AMANZI_FLOW_RELATIONS_SOIL_RESISTANCE_PARTITION_

#include "soil_resistance_sakagucki_zeng_model.hh"
#include "MeshPartition.hh"

namespace Amanzi {
namespace Flow {

typedef std::vector<Teuchos::RCP<SoilResistanceSakaguckiZengModel>>
  SoilResistanceSakaguckiZengModelList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, SoilResistanceSakaguckiZengModelList>
  SoilResistanceModelPartition;

// Non-member factory
Teuchos::RCP<SoilResistanceModelPartition>
createSoilResistanceModelPartition(Teuchos::ParameterList& plist);

} // namespace Flow
} // namespace Amanzi

#endif
