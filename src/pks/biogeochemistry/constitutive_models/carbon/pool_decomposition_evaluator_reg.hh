/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "pool_decomposition_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace BGC {
namespace BGCRelations {

// registry of method
Utils::RegisteredFactory<Evaluator, PoolDecompositionEvaluator> PoolDecompositionEvaluator::fac_(
  "pool decomposition evaluator");

} // namespace BGC
} // namespace ATS_Physics
} // namespace BGC
} // namespace ATS_Physics
} // namespace Amanzi
