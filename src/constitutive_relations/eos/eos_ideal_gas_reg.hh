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

#include "eos_ideal_gas.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

// registry of method
Utils::RegisteredFactory<EOS, EOSIdealGas> EOSIdealGas::factory_("ideal gas");

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
