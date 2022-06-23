/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  SubgridDisaggregateEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "SubgridDisaggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,SubgridDisaggregateEvaluator> SubgridDisaggregateEvaluator::factory_("subgrid disaggregate evaluator");

} // namespace
} // namespace
