#include "richards_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

Utils::RegisteredFactory<Evaluator, RichardsEnergyEvaluator>
  RichardsEnergyEvaluator::reg_("richards energy");

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
