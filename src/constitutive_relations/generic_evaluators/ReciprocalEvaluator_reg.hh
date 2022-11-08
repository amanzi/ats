/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ReciprocalEvaluator is the generic evaluator for dividing two vectors.

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/

#include "ReciprocalEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,ReciprocalEvaluator> ReciprocalEvaluator::factory_("reciprocal evaluator");

} // namespace
} // namespace
