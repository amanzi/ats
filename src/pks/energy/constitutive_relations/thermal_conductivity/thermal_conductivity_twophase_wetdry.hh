/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.

See ATS process model documentation's permafrost model for details.

Usage:

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="Thermal Conductivity Type" type="string" value="two-phase wet/dry"/>

    <Parameter name="thermal conductivity, wet" type="double" value=""/>
    <Parameter name="thermal conductivity, dry" type="double" value=""/>

    <Parameter name="epsilon" type="double" value="1.e-10"/>
    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
  </ParameterList>

Units: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_WETDRY_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_WETDRY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_twophase.hh"

namespace Amanzi {
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
} // namespace Amanzi

#endif
