/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Three-phase thermal conductivity based on paper by Peters-Lidard.
/*!

A three-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Empirical relationship for frozen soil based on Peters-Lidard

See Atchley et al GMD 2015 Supplementary Material for equations.

`"thermal conductivity type`" = `"three-phase wet/dry`"

.. _thermal_conductivity_threephase_wetdry-spec:
.. admonition:: thermal_conductivity_threephase_wetdry-spec

   * `"region`" ``[string]`` Region name on which to apply these parameters.
   * `"thermal conductivity, saturated (unfrozen) [W m^-1 K^-1]`" ``[double]``
     Thermal conductivity of fully saturated, unfrozen bulk soil.
   * `"thermal conductivity, dry [W m^-1 K^-1]`" ``[double]`` Thermal conductivity
     of fully dried bulk soil.
   * `"unsaturated alpha unfrozen [-]`" ``[double]`` Interpolating exponent
   * `"unsaturated alpha frozen [-]`" ``[double]`` Interpolating exponent
   * `"saturated beta frozen [-]`" ``[double]`` **1.0** Interpolating exponent
   * `"epsilon`" ``[double]`` **1e-10** Epsilon to keep saturations bounded away from 0.

*/

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_WETDRY_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_WETDRY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_threephase.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityThreePhaseWetDry : public ThermalConductivityThreePhase {
 public:
  ThermalConductivityThreePhaseWetDry(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq, double sat_ice, double temp);
  double
  DThermalConductivity_DPorosity(double porosity, double sat_liq, double sat_ice, double temp);
  double DThermalConductivity_DSaturationLiquid(double porosity,
                                                double sat_liq,
                                                double sat_ice,
                                                double temp);
  double
  DThermalConductivity_DSaturationIce(double porosity, double sat_liq, double sat_ice, double temp);
  double
  DThermalConductivity_DTemperature(double porosity, double sat_liq, double sat_ice, double temp);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_f_;
  double alpha_u_;
  double k_sat_u_;
  double k_dry_;
  double beta_sat_f_;

 private:
  static Utils::RegisteredFactory<ThermalConductivityThreePhase,
                                  ThermalConductivityThreePhaseWetDry>
    factory_;
};

} // namespace Energy
} // namespace Amanzi

#endif
