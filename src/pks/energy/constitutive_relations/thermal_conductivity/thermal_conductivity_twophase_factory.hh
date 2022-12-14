/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_TC_TWOPHASE_FACTORY_HH_
#define PK_ENERGY_RELATIONS_TC_TWOPHASE_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "thermal_conductivity_twophase.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityTwoPhaseFactory : public Utils::Factory<ThermalConductivityTwoPhase> {
 public:
  Teuchos::RCP<ThermalConductivityTwoPhase>
  createThermalConductivityModel(Teuchos::ParameterList& plist);
};

} // namespace Energy
} // namespace Amanzi

#endif
