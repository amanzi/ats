/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  The settlement evaluator gets the erosion rates.

*/

#include "settlement_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator, SettlementRateEvaluator>
  SettlementRateEvaluator ::factory_("settlement rate");

} // namespace Amanzi
