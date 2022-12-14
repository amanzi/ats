/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "wrm_linear_system.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<WRM, WRMLinearSystem> WRMLinearSystem::factory_("linear system");

} // namespace Flow
} // namespace Amanzi
