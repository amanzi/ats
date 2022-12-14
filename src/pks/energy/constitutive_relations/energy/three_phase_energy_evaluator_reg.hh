#include "three_phase_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

Utils::RegisteredFactory<Evaluator, ThreePhaseEnergyEvaluator>
  ThreePhaseEnergyEvaluator::reg_("three phase energy");

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
