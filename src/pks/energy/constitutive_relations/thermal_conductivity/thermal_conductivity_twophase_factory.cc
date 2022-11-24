/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for TC implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "thermal_conductivity_twophase_factory.hh"

namespace Amanzi {
namespace Energy {

// method for instantiating implementations
Teuchos::RCP<ThermalConductivityTwoPhase>
ThermalConductivityTwoPhaseFactory::createThermalConductivityModel(Teuchos::ParameterList& plist)
{
  std::string tc_typename = plist.get<std::string>("thermal conductivity type");
  return Teuchos::rcp(CreateInstance(tc_typename, plist));
};

} // namespace Energy
} // namespace Amanzi
