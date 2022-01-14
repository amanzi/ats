/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The settlement evaluator gets the erosion rates.

  License: BSD
  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "settlement_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator,SettlementRateEvaluator > SettlementRateEvaluator ::factory_("settlement rate");

} // namespace
