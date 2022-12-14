/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the flow component of the Amanzi code.
*/


#include "wrm_interfrost.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMInterfrost> WRMInterfrost::factory_("interfrost wrm");

} // namespace Flow
} // namespace Amanzi
