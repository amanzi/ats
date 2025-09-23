/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "wrm_van_genuchten.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMVanGenuchten> WRMVanGenuchten::factory_("van Genuchten");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
