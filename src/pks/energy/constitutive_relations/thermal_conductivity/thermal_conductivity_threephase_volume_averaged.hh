/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A volume-averaged thermal conductivity based on TCs of raw components.
/*!

A simple model of three-phase thermal conductivity, based upon volume-averaging
of four consitutive components.

`"thermal conductivity type`" = `"three-phase volume averaged`"

See Atchley et al GMD 2015 Supplementary Material for equations.

.. _thermal-conductivity-volume-averaged-spec:
.. admonition:: thermal-conductivity-volume-averaged-spec

    * `"thermal conductivity of soil [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of soil **grains**
    * `"thermal conductivity of liquid [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of liquid water.
    * `"thermal conductivity of gas [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of air.
    * `"thermal conductivity of ice [W m^-1 K^-1]`" ``[double]`` Thermal
      conductivity of frozen water.

*/

#pragma once

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_threephase.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityThreePhaseVolumeAveraged : public ThermalConductivityThreePhase {
 public:
  ThermalConductivityThreePhaseVolumeAveraged(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq, double sat_ice, double temp);
  double DThermalConductivity_DPorosity(double porosity, double sat_liq, double sat_ice, double temp);
  double DThermalConductivity_DSaturationLiquid(double porosity, double sat_liq, double sat_ice, double temp);
  double DThermalConductivity_DSaturationIce(double porosity, double sat_liq, double sat_ice, double temp);
  double DThermalConductivity_DTemperature(double porosity, double sat_liq, double sat_ice, double temp);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double k_soil_;
  double k_ice_;
  double k_liquid_;
  double k_gas_;

 private:
  static Utils::RegisteredFactory<ThermalConductivityThreePhase,
                                  ThermalConductivityThreePhaseVolumeAveraged>
    factory_;
};

} // namespace Energy
} // namespace Amanzi
