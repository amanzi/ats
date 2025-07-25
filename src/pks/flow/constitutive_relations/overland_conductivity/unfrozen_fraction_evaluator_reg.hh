/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the conductivity of surface flow.

*/

#include "unfrozen_fraction_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, UnfrozenFractionEvaluator> UnfrozenFractionEvaluator::fac_(
  "unfrozen fraction");

} // namespace Flow
} // namespace Amanzi
