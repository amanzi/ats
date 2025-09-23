/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "wrm_linear_system.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// registry of method
Utils::RegisteredFactory<WRM, WRMLinearSystem> WRMLinearSystem::factory_("linear system");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
