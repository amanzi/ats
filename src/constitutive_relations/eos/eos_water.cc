/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_water.hh"

namespace Amanzi {
namespace Relations {

EOSWater::EOSWater(Teuchos::ParameterList& eos_plist)
  : EOSConstantMolarMass(0.0180153),
    eos_plist_(eos_plist),
    ka_(999.915),
    kb_(0.0416516),
    kc_(-0.0100836),
    kd_(0.000206355),
    kT0_(273.15),
    kalpha_(5.0e-10),
    kp0_(1.0e5)
{
  if (eos_plist_.isParameter("molar mass [kg mol^-1]")) {
    M_ = eos_plist_.get<double>("molar mass [kg mol^-1]");
  } else {
    M_ = eos_plist_.get<double>("molar mass [g mol^-1]", 18.0153) * 1e-3;
  }
};


double
EOSWater::MassDensity(std::vector<double>& params)
{
  double T = params[0];
  double p = std::max(params[1], 101325.);

  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_ * dT) * dT) * dT;
  return rho1bar * (1.0 + kalpha_ * (p - kp0_));
};


double
EOSWater::DMassDensityDT(std::vector<double>& params)
{
  double T = params[0];
  double p = std::max(params[1], 101325.);

  double dT = T - kT0_;
  double rho1bar = kb_ + (2.0 * kc_ + 3.0 * kd_ * dT) * dT;
  return rho1bar * (1.0 + kalpha_ * (p - kp0_));
};


double
EOSWater::DMassDensityDp(std::vector<double>& params)
{
  double T = params[0];
  double p = params[1];

  if (p < 101325.) {
    return 0.;
  } else {
    double dT = T - kT0_;
    double rho1bar = ka_ + (kb_ + (kc_ + kd_ * dT) * dT) * dT;
    return rho1bar * kalpha_;
  }
};

} // namespace Relations
} // namespace Amanzi
