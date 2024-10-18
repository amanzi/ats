/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

*/

#include "errors.hh"
#include "viscosity_water.hh"

namespace Amanzi {
namespace Relations {

ViscosityWater::ViscosityWater(Teuchos::ParameterList& eos_plist)
  : eos_plist_(eos_plist),
    kav1_(998.333),
    kbv1_(-8.1855),
    kcv1_(0.00585),
    kbv2_(1.3272),
    kcv2_(-0.001053),
    kT1_(293.15){};


double
ViscosityWater::Viscosity(double T)
{
  double dT = kT1_ - T;
  double xi;
  if (T < kT1_) {
    double A = kav1_ + (kbv1_ + kcv1_ * dT) * dT;
    xi = 1301.0 * (1.0 / A - 1.0 / kav1_);
  } else {
    double A = (kbv2_ + kcv2_ * dT) * dT;
    xi = A / (T - 168.15);
  }
  double visc = 0.001 * std::pow(10.0, xi);

  if (visc < 1.e-16) {
    std::cout << "Invalid temperature, T = " << T << std::endl;
    Exceptions::amanzi_throw(Errors::CutTimestep());
  }
  return visc;
};

double
ViscosityWater::DViscosityDT(double T)
{
  double dT = kT1_ - T;
  double xi;
  double dxi_dT;

  if (T < kT1_) {
    double A = kav1_ + (kbv1_ + kcv1_ * dT) * dT;
    double dxi_dA = -1301. * std::pow(A, -2);
    double dA_dT = -(kbv1_ + 2 * kcv1_ * dT);
    xi = 1301.0 * (1.0 / A - 1.0 / kav1_);
    dxi_dT = dxi_dA * dA_dT;

  } else {
    double A = (kbv2_ + kcv2_ * dT) * dT;
    double dA_dT = kbv2_ + 2 * kcv2_ * dT;
    xi = A / (T - 168.15);
    dxi_dT = dA_dT / (T - 168.15) - A * std::pow(T - 168.15, -2);
  }


  double dvisc_dxi = 0.001 * std::pow(10.0, xi) * std::log(10.);
  return dvisc_dxi * dxi_dT;
};


} // namespace Relations
} // namespace Amanzi
