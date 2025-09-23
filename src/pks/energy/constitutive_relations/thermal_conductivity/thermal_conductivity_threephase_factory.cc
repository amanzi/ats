/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for TC implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "thermal_conductivity_threephase_factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

// method for instantiating implementations
Teuchos::RCP<ThermalConductivityThreePhase>
ThermalConductivityThreePhaseFactory::createThermalConductivityModel(Teuchos::ParameterList& plist)
{
  std::string tc_typename = plist.get<std::string>("thermal conductivity type");
  return Teuchos::rcp(CreateInstance(tc_typename, plist));
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
