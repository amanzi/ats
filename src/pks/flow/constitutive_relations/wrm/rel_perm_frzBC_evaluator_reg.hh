/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov), Bo Gao (gaob@ornl.gov)
*/

#include "rel_perm_frzBC_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, RelPermFrzBCEvaluator>
  RelPermFrzBCEvaluator::factory_("relative permeability, freezing Brooks-Corey");

} // namespace Flow
} // namespace Amanzi
