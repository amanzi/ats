/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  PCLiquidEvaluator is the interface between state/data and the model, a PC relation.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "pc_liq_atm.hh"
#include "pc_liquid_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, PCLiquidEvaluator>
  PCLiquidEvaluator::factory_("capillary pressure, atmospheric gas over liquid");

} // namespace Flow
} // namespace Amanzi
