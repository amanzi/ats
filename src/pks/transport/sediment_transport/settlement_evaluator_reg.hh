/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

/*
  The settlement evaluator gets the erosion rates.

  License: BSD
  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "settlement_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator, SettlementRateEvaluator>
  SettlementRateEvaluator ::factory_("settlement rate");

} // namespace Amanzi
