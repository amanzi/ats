/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, MultiplicativeEvaluator>
  MultiplicativeEvaluator::factory_("multiplicative evaluator");

} // namespace Relations
} // namespace Amanzi
