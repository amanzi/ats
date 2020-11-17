/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include <cmath>
#include "thermal_conductivity_threephase_peterslidard.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivityThreePhasePetersLidard::ThermalConductivityThreePhasePetersLidard(
      Teuchos::ParameterList& plist) : plist_(plist) {
  InitializeFromPlist_();
};

double ThermalConductivityThreePhasePetersLidard::ThermalConductivity(double poro,
        double sat_liq, double sat_ice, double temp) {
  double k_dry = (d_*(1-poro)*k_soil_ + k_gas_*poro)/(d_*(1-poro) + poro);
  double k_sat_u = pow(k_soil_,(1-poro)) * pow(k_liquid_,poro);
  double k_sat_f = pow(k_soil_,(1-poro)) * pow(k_ice_,poro);
  double kersten_u = pow(sat_liq + eps_, alpha_u_);
  double kersten_f = pow(sat_ice + eps_, alpha_f_);
  return kersten_f * k_sat_f + kersten_u * k_sat_u
    + (1.0 - kersten_f - kersten_u) * k_dry;
};

void ThermalConductivityThreePhasePetersLidard::InitializeFromPlist_() {
  d_ = 0.053; // unitless empericial parameter

  eps_ = plist_.get<double>("epsilon [-]", 1.e-10);
  alpha_u_ = plist_.get<double>("unsaturated alpha unfrozen [-]");
  alpha_f_ = plist_.get<double>("unsaturated alpha frozen [-]");
  k_soil_ = plist_.get<double>("thermal conductivity of soil [W m^-1 K^-1]");
  k_ice_ = plist_.get<double>("thermal conductivity of ice [W m^-1 K^-1]");
  k_liquid_ = plist_.get<double>("thermal conductivity of liquid [W m^-1 K^-1]");
  k_gas_ = plist_.get<double>("thermal conductivity of gas [W m^-1 K^-1]");
};

} // namespace Relations
} // namespace Energy
