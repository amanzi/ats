/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  PCLiquidEvaluator is the interface between state/data and the model, a PC relation.

*/

#include "pc_liq_atm.hh"
#include "pc_liquid_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, PCLiquidEvaluator> PCLiquidEvaluator::factory_(
  "capillary pressure, atmospheric gas over liquid");

} // namespace Flow
} // namespace Amanzi
