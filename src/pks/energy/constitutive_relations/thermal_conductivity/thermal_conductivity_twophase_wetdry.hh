/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.

`"thermal conductivity type`" = `"two-phase wet/dry`"

.. _thermal-conductivity-twophase-wetdry-spec:
.. admonition:: thermal-conductivity-twophase-wetdry-spec

   * `"thermal conductivity, wet [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of saturated soil
   * `"thermal conductivity, dry [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of dry soil
   * `"unsaturated alpha [-]`" ``[double]`` Interpolating exponent
   * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

Example:

.. code:: xml

  <ParameterList name="thermal conductivity model">
    <Parameter name="thermal conductivity type" type="string" value="two-phase wet/dry"/>
    <Parameter name="thermal conductivity, wet [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity, dry [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="epsilon" type="double" value="1.e-10"/>
    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
  </ParameterList>

Units: [W m^-1 K^-1]

*/

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_WETDRY_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_WETDRY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_twophase.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class ThermalConductivityTwoPhaseWetDry : public ThermalConductivityTwoPhase {
 public:
  ThermalConductivityTwoPhaseWetDry(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_;
  double k_wet_;
  double k_dry_;

 private:
  static Utils::RegisteredFactory<ThermalConductivityTwoPhase, ThermalConductivityTwoPhaseWetDry>
    factory_;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
