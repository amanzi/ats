/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "rel_perm_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, RelPermEvaluator> RelPermEvaluator::factory_(
  "relative permeability, water retention model");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
