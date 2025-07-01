/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  The trapping evaluator gets the erosion rates.

*/

#include "trapping_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator, TrappingRateEvaluator>
  TrappingRateEvaluator ::factory_("trapping rate");

} // namespace Amanzi
