/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "surface_relperm_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, SurfaceRelPermEvaluator>
  SurfaceRelPermEvaluator::fac_("surface rel perm");

} // namespace Flow
} // namespace Amanzi
