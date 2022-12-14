/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

  Authors: Ethan Coon (ecoon@lanl.gov)
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
  R_ = eos_plist_.get<double>("Ideal gas constant [J/mol-K]", 8.3144621);

  if (eos_plist_.isParameter("Molar mass of gas [kg/mol]")) {
    M_ = eos_plist_.get<double>("Molar mass of gas [kg/mol]");
  } else {
    M_ = eos_plist_.get<double>("Molar mass of gas [g/mol]", 28.956) * 1e-3;
  }
};

} // namespace Relations
} // namespace Amanzi
