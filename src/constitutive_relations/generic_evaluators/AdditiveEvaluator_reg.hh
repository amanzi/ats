/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  AdditiveEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "AdditiveEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,AdditiveEvaluator> AdditiveEvaluator::factory_("additive evaluator");

} // namespace
} // namespace
