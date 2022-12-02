/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

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
