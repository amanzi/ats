/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "wrm_interfrost.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMInterfrost> WRMInterfrost::factory_("interfrost wrm");

} // namespace Flow
} // namespace Amanzi
