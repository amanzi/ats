/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Internal energy based on a quadratic fit to data.

/*!

Quadratic internal energy model -- function of Cv and temperature

.. math::

    u = L_f + C_v * (T - T_{ref}) + b(T - T_{ref})^2

*/
#include "iem_quadratic.hh"

namespace Amanzi {
namespace Energy {

IEMQuadratic::IEMQuadratic(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double IEMQuadratic::InternalEnergy(double temp) {
  double dT = temp - T0_;
  return u0_ + (ka_ + kb_*dT) * dT;
};

double IEMQuadratic::DInternalEnergyDT(double temp) {
  double dT = temp - T0_;
  return ka_ + 2.0*kb_*dT;
};

void IEMQuadratic::InitializeFromPlist_() {
  if (plist_.isParameter("quadratic u_0 [J kg^-1]")) {
    u0_ = 1.e-6 * plist_.get<double>("latent heat [J kg^-1]");
    ka_ = 1.e-6 * plist_.get<double>("heat capacity [J kg^-1 K^-1]");
    kb_ = 1.e-6 * plist_.get<double>("quadratic b [J kg^-1 K^-2]");
    molar_basis_ = false;

  } else {
    u0_ = 1.e-6 * plist_.get<double>("latent heat [J mol^-1]");
    ka_ = 1.e-6 * plist_.get<double>("heat capacity [J mol^-1 K^-1]");
    kb_ = 1.e-6 * plist_.get<double>("quadratic b [J mol^-1 K^-2]");
    molar_basis_ = true;
  }

  T0_ = plist_.get<double>("reference temperature [K]", 273.15);
};

} // namespace
} // namespace
