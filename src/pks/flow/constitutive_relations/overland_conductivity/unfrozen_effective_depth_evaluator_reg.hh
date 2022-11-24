/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen effective depth.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_effective_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, UnfrozenEffectiveDepthEvaluator>
  UnfrozenEffectiveDepthEvaluator::fac_("unfrozen effective depth");

} // namespace Flow
} // namespace Amanzi
