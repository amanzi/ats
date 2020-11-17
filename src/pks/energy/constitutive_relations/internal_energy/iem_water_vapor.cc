/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Internal energy model for air and water vapor.

/*!

Internal energy model for air and water vapor, relative to water @237.15K

*/
#include "iem_water_vapor.hh"

namespace Amanzi {
namespace Energy {

IEMWaterVapor::IEMWaterVapor(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double IEMWaterVapor::InternalEnergy(double temp, double mol_frac_gas) {
  return (1.0 + 0.622*mol_frac_gas) * Cv_air_ * (temp - 273.15) + mol_frac_gas*heat_vaporization_;
};

double IEMWaterVapor::DInternalEnergyDT(double temp, double mol_frac_gas) {
  // evaluated at constant mol_frac gas for now?
  return (1.0 + 0.622*mol_frac_gas) * Cv_air_;
};

double IEMWaterVapor::DInternalEnergyDomega(double temp, double mol_frac_gas) {
  // evaluated at constant mol_frac gas for now?
  return heat_vaporization_ + 0.622 * Cv_air_ * (temp - 273.15);
};

void IEMWaterVapor::InitializeFromPlist_() {
  Cv_air_ = 1.e-6 * plist_.get<double>("heat capacity [J mol^-1 K^-1]", 13.0);
  heat_vaporization_ = 1.e-6 * plist_.get<double>("latent heat [J mol^-1]", 4.065e4);
};

} // namespace
} // namespace
