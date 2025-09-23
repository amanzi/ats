/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Three-phase thermal conductivity models are used to compute the bulk thermal
conductivity from the constitutive parts -- gas, liquid, ice, and background
material.

.. _thermal-conductivity-threephase-typed-spec:
.. admonition:: thermal-conductivity-threephase-typed-spec

   * `"region`" ``[string]`` Region on which the model is valid.
   * `"thermal conductivity type`" ``[string]`` Name of the model, see below for options.
   * `"_thermal_conductivity_type_ parameters`"
     ``[_thermal_conductivity_type_-spec]`` See below for the required
     parameter spec for each type.

*/

#ifndef PK_ENERGY_RELATIONS_TC_THREEPHASE_FACTORY_HH_
#define PK_ENERGY_RELATIONS_TC_THREEPHASE_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "thermal_conductivity_threephase.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class ThermalConductivityThreePhaseFactory : public Utils::Factory<ThermalConductivityThreePhase> {
 public:
  Teuchos::RCP<ThermalConductivityThreePhase> createThermalConductivityModel(
    Teuchos::ParameterList& plist);
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
