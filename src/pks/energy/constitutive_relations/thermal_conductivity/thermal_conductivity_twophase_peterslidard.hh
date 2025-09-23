/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

A two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Emperical fit for dry conductivity from Peters-Lidard et al '98.

See Atchley et al GMD 2015 Supplementary Material for equations.

`"thermal conductivity type`" = `"two-phase Peters-Lidard`"

.. _thermal-conductivity-twophase-peterslidard-spec:
.. admonition:: thermal-conductivity-twophase-peterslidard-spec

    * `"region`" ``[string]`` Region name on which to apply these parameters.
    * `"thermal conductivity of soil [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of soil grains (not bulk soil)
    * `"thermal conductivity of liquid [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of liquid (water)
    * `"thermal conductivity of gas [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of gas (air)
    * `"unsaturated alpha [-]`" ``[double]`` Interpolating exponent
    * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

Example:

.. code:: xml

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="thermal conductivity type" type="string" value="two-phase Peters-Lidard"/>
    <Parameter name="thermal conductivity of soil [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of gas [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
    <Parameter name="epsilon" type="double" value="1.e-10"/>
  </ParameterList>

Units: [W m^-1 K^-1]

*/
#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_PETERSLIDARD_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_PETERSLIDARD_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_twophase.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class ThermalConductivityTwoPhasePetersLidard : public ThermalConductivityTwoPhase {
 public:
  ThermalConductivityTwoPhasePetersLidard(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_;
  double k_soil_;
  double k_liquid_;
  double k_gas_;
  double d_;

 private:
  static Utils::RegisteredFactory<ThermalConductivityTwoPhase,
                                  ThermalConductivityTwoPhasePetersLidard>
    factory_;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
