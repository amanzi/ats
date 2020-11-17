/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include <cmath>
#include "thermal_conductivity_threephase_sutra_hacked.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivityThreePhaseSutraHacked::ThermalConductivityThreePhaseSutraHacked(
      Teuchos::ParameterList& plist) : plist_(plist) {
  InitializeFromPlist_();
};

double ThermalConductivityThreePhaseSutraHacked::ThermalConductivity(double poro,
        double sat_liq, double sat_ice, double temp) {
  if (sat_ice == 0.) {
    return k_unfrozen_;
  } else if (std::abs(sat_ice - (1-sr_)) < 1.e-10) {
    return k_frozen_;
  }
  return k_mushy_;
};

void ThermalConductivityThreePhaseSutraHacked::InitializeFromPlist_() {
  k_frozen_ = plist_.get<double>("thermal conductivity of frozen zone [W m^-1 K^-1]");
  k_unfrozen_ = plist_.get<double>("thermal conductivity of unfrozen zone [W m^-1 K^-1]");
  k_mushy_ = plist_.get<double>("thermal conductivity of mushy zone [W m^-1 K^-1]");
  sr_ = plist_.get<double>("residual saturation [-]");
};

} // namespace Relations
} // namespace Energy
