/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "rel_perm_sutraice_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, RelPermSutraIceEvaluator> RelPermSutraIceEvaluator::factory_(
  "relative permeability, SutraICE");

} // namespace Flow
} // namespace Amanzi
