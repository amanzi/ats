/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Three-phase thermal conductivity based on paper by Peters-Lidard.

/*!

A three-phase thermal conductivity, based upon:

- A mixture model using interpolation across various components.
- Power-law Kersten number.

See Atchley et al GMD 2015 Supplementary Material for equations.

.. _thermal-conductivity-threephase-peterslidard-spec:
.. admonition:: thermal-conductivity-threephase-peterslidard-spec

    * `"thermal conductivity of soil [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of soil grains (not bulk soil)
    * `"thermal conductivity of liquid [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of liquid (water)
    * `"thermal conductivity of gas [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of gas (air)
    * `"thermal conductivity of ice [W m^-1 K^-1]`" ``[double]`` Thermal conductivity of ice
    * `"unsaturated alpha unfrozen [-]`" ``[double]`` Interpolating exponent
    * `"unsaturated alpha frozen [-]`" ``[double]`` Interpolating exponent
    * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

Usage:

  <ParameterList name="thermal_conductivity">
    <Parameter name="thermal conductivity type" type="string" value="three-phase Peters-Lidard"/>
    <Parameter name="thermal conductivity of soil [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of gas [W m^-1 K^-1]" type="double" value=""/>
    <Parameter name="thermal conductivity of ice [W m^-1 K^-1]" type="double" value=""/>

    <Parameter name="unsaturated alpha unfrozen [-]" type="double" value=""/>
    <Parameter name="unsaturated alpha frozen [-]" type="double" value=""/>

    <Parameter name="epsilon" type="double" value="1.e-10"/>
  </ParameterList>

*/

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_PETERSLIDARD_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_PETERSLIDARD_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_threephase.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityThreePhasePetersLidard : public ThermalConductivityThreePhase {

public:
  ThermalConductivityThreePhasePetersLidard(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq, double sat_ice, double temp);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_f_;
  double alpha_u_;
  double k_soil_;
  double k_ice_;
  double k_liquid_;
  double k_gas_;
  double d_;

private:
  static Utils::RegisteredFactory<ThermalConductivityThreePhase,
                                  ThermalConductivityThreePhasePetersLidard> factory_;

};

}
}

#endif
