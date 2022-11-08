/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  InitialTimeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "InitialTimeEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,InitialTimeEvaluator> InitialTimeEvaluator::reg_("initial value");

} // namespace
} // namespace
