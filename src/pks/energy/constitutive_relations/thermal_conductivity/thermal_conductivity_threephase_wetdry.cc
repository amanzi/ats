/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include <cmath>
#include "thermal_conductivity_threephase_wetdry.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivityThreePhaseWetDry::ThermalConductivityThreePhaseWetDry(
      Teuchos::ParameterList& plist) : plist_(plist) {
  InitializeFromPlist_();
};

double ThermalConductivityThreePhaseWetDry::ThermalConductivity(double poro,
        double sat_liq, double sat_ice, double temp) {
  double Ki = 831.51 * std::pow(temp, -1.0552);
  double Kl = 0.5611;
  double k_sat_f =  beta_sat_f_ * k_sat_u_ * std::pow(Ki/Kl, poro);

  double kersten_u = std::pow(sat_liq + eps_, alpha_u_);
  double kersten_f = std::pow(sat_ice + eps_, alpha_f_);
  return kersten_f * k_sat_f + kersten_u * k_sat_u_ + (1.0 - kersten_f - kersten_u) * k_dry_;
};

double
ThermalConductivityThreePhaseWetDry::DThermalConductivity_DPorosity(double poro,
        double sat_liq, double sat_ice, double temp) {
  double Ki = 831.51 * std::pow(temp, -1.0552);
  double Kl = 0.5611;
  double kersten_f = std::pow(sat_ice + eps_, alpha_f_);
  double dk_sat_f =  beta_sat_f_ * k_sat_u_ * std::pow(Ki/Kl, poro) * std::log(Ki/Kl);
  return kersten_f * dk_sat_f;
}

double
ThermalConductivityThreePhaseWetDry::DThermalConductivity_DSaturationLiquid(double poro,
        double sat_liq, double sat_ice, double temp) {
  double dkersten_u = alpha_u_ * std::pow(sat_liq + eps_, alpha_u_ - 1.0);
  return dkersten_u * k_sat_u_ - dkersten_u * k_dry_;

}

double
ThermalConductivityThreePhaseWetDry::DThermalConductivity_DSaturationIce(double poro,
        double sat_liq, double sat_ice, double temp) {
  double Ki = 831.51 * std::pow(temp, -1.0552);
  double Kl = 0.5611;
  double k_sat_f =  beta_sat_f_ * k_sat_u_ * std::pow(Ki/Kl, poro);
  double dkersten_f = alpha_f_ * std::pow(sat_ice + eps_, alpha_f_-1.0);
  return dkersten_f * k_sat_f - dkersten_f * k_dry_;
}

double
ThermalConductivityThreePhaseWetDry::DThermalConductivity_DTemperature(double poro,
        double sat_liq, double sat_ice, double temp) {
  double Ki = 831.51 * std::pow(temp, -1.0552);
  double Kl = 0.5611;
  double dKi = 831.51 * -1.0552 * std::pow(temp, -2.0552);
  double dk_sat_f =  beta_sat_f_ * k_sat_u_ * poro * std::pow(Ki/Kl, poro-1.0) * dKi / Kl;
  double kersten_f = std::pow(sat_ice + eps_, alpha_f_);
  return kersten_f * dk_sat_f;
}


void ThermalConductivityThreePhaseWetDry::InitializeFromPlist_() {
  eps_ = plist_.get<double>("epsilon [-]", 1.e-10);
  alpha_u_ = plist_.get<double>("unsaturated alpha unfrozen [-]");
  alpha_f_ = plist_.get<double>("unsaturated alpha frozen [-]");
  k_dry_ = plist_.get<double>("thermal conductivity, dry [W m^-1 K^-1]");
  k_sat_u_ = plist_.get<double>("thermal conductivity, saturated (unfrozen) [W m^-1 K^-1]");
  beta_sat_f_ = plist_.get<double>("saturated beta frozen [-]",1.0);
};

} // namespace Relations
} // namespace Energy
