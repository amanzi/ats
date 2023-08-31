/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the unfrozen effective depth.

*/

#include "unfrozen_effective_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, UnfrozenEffectiveDepthEvaluator>
  UnfrozenEffectiveDepthEvaluator::fac_("unfrozen effective depth");

} // namespace Flow
} // namespace Amanzi
