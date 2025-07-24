/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  ATS

  EOS for salt water (does not implement viscosity at this point!)
  For this model dependcy on concentration is only assumed.

*/

#include "eos_factory.hh"
#include "eos_sw.hh"


namespace Amanzi {
namespace Relations {

EOS_SW::EOS_SW(Teuchos::ParameterList& eos_plist)
  : eos_plist_(eos_plist)
{
  InitializeFromPlist_();
}

double
EOS_SW::MassDensity(std::vector<double>& params)
{
  double C = params[0];
  return rho_f_ + E_ * C;
};

double
EOS_SW::DMassDensityDMoleFraction(std::vector<double>& params)
{
  return E_;
};

double
EOS_SW::MolarDensity(std::vector<double>& params)
{
  double C = params[0];
  return MassDensity(params) / (M_water_ * (1 - C) + M_salt_ * C);
};

double
EOS_SW::DMolarDensityDMoleFraction(std::vector<double>& params)
{
  double C = params[0];
  double b = (M_water_ * (1 - C) + M_salt_ * C);
  return (DMassDensityDMoleFraction(params) * b - MassDensity(params) * (M_salt_ - M_water_)) /
         (b * b);
};


void
EOS_SW::InitializeFromPlist_()
{
  E_ = eos_plist_.get<double>("sea water density coefficient", 750);
  rho_f_ = eos_plist_.get<double>("fresh water mass density [kg/m^3]", 1000.0);
  if (eos_plist_.isParameter("salt molar mass [kg/mol]")) {
    M_salt_ = eos_plist_.get<double>("salt molar mass [kg/mol]");
  } else {
    M_salt_ = eos_plist_.get<double>("salt molar mass [g/mol]", 58.5) * 1e-3;
  }

  if (eos_plist_.isParameter("water molar mass [kg/mol]")) {
    M_water_ = eos_plist_.get<double>("water molar mass [kg/mol]");
  } else {
    M_water_ = eos_plist_.get<double>("water molar mass [g/mol]", 18.0153) * 1e-3;
  }
};

} // namespace Relations
} // namespace Amanzi
