/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  EOS for liquid ice.

*/

#include "eos_ice.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

EOSIce::EOSIce(Teuchos::ParameterList& eos_plist)
  : eos_plist_(eos_plist),
    ka_(916.724),
    kb_(-0.147143),
    kc_(-0.000238095),
    kT0_(273.15),
    kalpha_(1.0e-10),
    kp0_(1.0e5)
{
  InitializeFromPlist_();
};

double
EOSIce::MassDensity(std::vector<double>& params)
{
  double T = params[0];
  double p = std::max(params[1], 101325.);

  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + kc_ * dT) * dT;
  return rho1bar * (1.0 + kalpha_ * (p - kp0_));
};


double
EOSIce::DMassDensityDT(std::vector<double>& params)
{
  double T = params[0];
  double p = std::max(params[1], 101325.);

  double dT = T - kT0_;
  double rho1bar = kb_ + 2.0 * kc_ * dT;
  return rho1bar * (1.0 + kalpha_ * (p - kp0_));
};


double
EOSIce::DMassDensityDp(std::vector<double>& params)
{
  double T = params[0];
  double p = params[1];

  if (p < 101325.) {
    return 0.;
  } else {
    double T = params[0];
    double p = params[1];
    double dT = T - kT0_;
    double rho1bar = ka_ + (kb_ + kc_ * dT) * dT;
    return rho1bar * kalpha_;
  }
};


void
EOSIce::InitializeFromPlist_()
{
  if (eos_plist_.isParameter("molar mass [kg mol^-1]")) {
    M_ = eos_plist_.get<double>("molar mass [kg mol^-1]");
  } else {
    M_ = eos_plist_.get<double>("molar mass [g mol^-1]", 18.0153) * 1e-3;
  }
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
