/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

#include "slope_stability_trigrs_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, SlopeStabilityTrigrsEvaluator>
  SlopeStabilityTrigrsEvaluator::fac_("slope stability, TRIGRS");

} // namespace Flow
} // namespace Amanzi
