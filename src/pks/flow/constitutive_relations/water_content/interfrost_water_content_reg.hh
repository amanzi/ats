#include "interfrost_water_content.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, InterfrostWaterContent>
  InterfrostWaterContent::reg_("interfrost water content");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
