/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.

  License: BSD
  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "erosion_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator, ErosionRateEvaluator>
  ErosionRateEvaluator ::factory_("erosion rate");

} // namespace Amanzi
