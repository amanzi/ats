/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "bioturbation_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

// registry of method
Utils::RegisteredFactory<Evaluator, BioturbationEvaluator> BioturbationEvaluator::fac_(
  "bioturbation evaluator");

} // namespace BGCRelations
} // namespace BGC
} // namespace Amanzi
