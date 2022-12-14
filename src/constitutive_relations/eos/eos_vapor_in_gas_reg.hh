/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_vapor_in_gas.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<EOS, EOSVaporInGas> EOSVaporInGas::factory_("vapor in gas");

} // namespace Relations
} // namespace Amanzi
