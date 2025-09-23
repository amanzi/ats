/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

*/

#include "eos_factory.hh"
#include "eos_vapor_in_gas.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

Utils::RegisteredFactory<EOS, EOSVaporInGas> EOSVaporInGas::factory_("vapor in gas");

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
