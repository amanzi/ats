/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  SubgridAggregateEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "SubgridAggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SubgridAggregateEvaluator> SubgridAggregateEvaluator::factory_("subgrid aggregate evaluator");

} // namespace
} // namespace
