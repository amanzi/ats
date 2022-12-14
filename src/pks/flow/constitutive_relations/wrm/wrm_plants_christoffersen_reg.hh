/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "wrm_plants_christoffersen.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMPlantChristoffersen>
  WRMPlantChristoffersen::factory_("plant Christoffersen");

} // namespace Flow
} // namespace Amanzi
