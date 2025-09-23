/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "subgrid_mobile_depth_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<Evaluator, SubgridMobileDepthEvaluator>
  SubgridMobileDepthEvaluator::factory_("subgrid mobile depth");

} // namespace FlowRelations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
