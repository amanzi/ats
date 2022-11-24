/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

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
