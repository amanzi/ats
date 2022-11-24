/*
  This is the flow component of the Amanzi code.
  License: BSD
*/


#include "wrm_constant.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMConstant> WRMConstant::factory_("constant");

} // namespace Flow
} // namespace Amanzi
