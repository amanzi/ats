/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

*/

#include "eos_ideal_gas.hh"

namespace Amanzi {
namespace Relations {

EOSIdealGas::EOSIdealGas(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist)
{
  InitializeFromPlist_();
};


double
EOSIdealGas::MolarDensity(std::vector<double>& params)
{
  double T = params[0];
  double p = std::max(params[1], 101325.);
  return p / (R_ * T);
};

double
EOSIdealGas::DMolarDensityDT(std::vector<double>& params)
{
  double T = params[0];
  double p = std::max(params[1], 101325.);
  return -p / (R_ * T * T);
};

double
EOSIdealGas::DMolarDensityDp(std::vector<double>& params)
{
  double T = params[0];
  double p = std::max(params[1], 101325.);
  return 1.0 / (R_ * T);
};


void
EOSIdealGas::InitializeFromPlist_()
{
  R_ = eos_plist_.get<double>("ideal gas constant [J mol^-1 K^-1]", 8.3144621);

  if (eos_plist_.isParameter("molar mass of gas [kg mol^-1]")) {
    M_ = eos_plist_.get<double>("molar mass of gas [kg mol^-1]");
  } else {
    M_ = eos_plist_.get<double>("molar mass of gas [g mol^-1]", 28.956) * 1e-3;
  }
};

} // namespace Relations
} // namespace Amanzi
