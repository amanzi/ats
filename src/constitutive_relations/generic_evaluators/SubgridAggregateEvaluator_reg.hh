/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

/*
  SubgridAggregateEvaluator is the generic evaluator for multipying two vectors.

*/

#include "SubgridAggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, SubgridAggregateEvaluator> SubgridAggregateEvaluator::factory_(
  "subgrid aggregate evaluator");

} // namespace Relations
} // namespace Amanzi
