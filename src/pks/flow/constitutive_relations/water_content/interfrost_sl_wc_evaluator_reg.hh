#include "interfrost_sl_wc_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, InterfrostSlWcEvaluator>
  InterfrostSlWcEvaluator::reg_("interfrost sl water content");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
