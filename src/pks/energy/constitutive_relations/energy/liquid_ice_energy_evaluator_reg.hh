#include "liquid_ice_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

Utils::RegisteredFactory<Evaluator, LiquidIceEnergyEvaluator>
  LiquidIceEnergyEvaluator::reg_("liquid+ice energy");

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
