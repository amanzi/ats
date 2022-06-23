/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.

  License: BSD
  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "organic_matter_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator,OrganicMatterRateEvaluator > OrganicMatterRateEvaluator ::factory_("organic matter rate");

} // namespace
