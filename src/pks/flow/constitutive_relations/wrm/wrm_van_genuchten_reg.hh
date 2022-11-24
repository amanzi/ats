/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "wrm_van_genuchten.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMVanGenuchten> WRMVanGenuchten::factory_("van Genuchten");

} // namespace Flow
} // namespace Amanzi
