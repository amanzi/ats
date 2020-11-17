/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Internal energy based on a linear fit.

/*!

Linear internal energy model -- function of Cv and temperature

.. math::

    u = L_f +  C_v * (T - T_{ref})
*/


#include "iem_linear.hh"

namespace Amanzi {
namespace Energy {

IEMLinear::IEMLinear(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double IEMLinear::InternalEnergy(double temp) {
  return L_ + Cv_ * (temp - T_ref_);
};

void IEMLinear::InitializeFromPlist_() {
  if (plist_.isParameter("heat capacity [J kg^-1 K^-1]")) {
    Cv_ = 1.e-6 * plist_.get<double>("heat capacity [J kg^-1 K^-1]");
    molar_basis_ = false;
  } else if (plist_.isParameter("heat capacity [MJ kg^-1 K^-1]")) {
    Cv_ = plist_.get<double>("heat capacity [MJ kg^-1 K^-1]");
    molar_basis_ = false;
  } else if (plist_.isParameter("heat capacity [MJ mol^-1 K^-1]")) {
    Cv_ = plist_.get<double>("heat capacity [MJ mol^-1 K^-1]");
    molar_basis_ = true;
  } else {
    Cv_ = 1.e-6 * plist_.get<double>("heat capacity [J mol^-1 K^-1]");
    molar_basis_ = true;
  }

  if (molar_basis_) {
    if (plist_.isParameter("latent heat [MJ mol^-1]")) {
      L_ = plist_.get<double>("latent heat [MJ mol^-1]");
    } else {
      L_ = 1.e-6 * plist_.get<double>("latent heat [J mol^-1]", 0.);
    }
  } else {
    if (plist_.isParameter("latent heat [MJ kg^-1]")) {
      L_ = plist_.get<double>("latent heat [MJ kg^-1]");
    } else {
      L_ = 1.e-6 * plist_.get<double>("latent heat [J kg^-1]", 0.);
    }
  }

  T_ref_ = plist_.get<double>("reference temperature [K]", 273.15);
};

} // namespace
} // namespace
