/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The trapping evaluator gets the erosion rates.

  License: BSD
  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "trapping_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator, TrappingRateEvaluator>
  TrappingRateEvaluator ::factory_("trapping rate");

} // namespace Amanzi
