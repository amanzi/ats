#include "interfrost_denergy_dtemperature_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, InterfrostDenergyDtemperatureEvaluator>
  InterfrostDenergyDtemperatureEvaluator::reg_("interfrost denergy_dtemperature");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
