/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "overland_conductivity_subgrid_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, OverlandConductivitySubgridEvaluator>
  OverlandConductivitySubgridEvaluator::factory_("overland conductivity subgrid");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
