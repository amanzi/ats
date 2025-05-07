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

EOS_SW::EOS_SW(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist)
{
  InitializeFromPlist_();
}

double
EOS_SW::MassDensity(std::vector<double>& params) // liquid
{
  // [kg liquid / m^3]
  double mol_ratio = params[0];
  return rho_f_ + E_ * M_salt_ * mol_ratio * MolarDensity(params);
};

double
EOS_SW::DMassDensityDMolarRatio(std::vector<double>& params)
{
  AMANZI_ASSERT(0);
};

double
EOS_SW::MolarDensity(std::vector<double>& params)
{
  // [mol H2O / m^3]
  double mol_ratio = params[0];
  return rho_f_ / (M_water_ + (1 - E_) * mol_ratio * M_salt_);
};

double
EOS_SW::DMolarDensityDMolarRatio(std::vector<double>& params)
{
  AMANZI_ASSERT(0);
};


void
EOS_SW::InitializeFromPlist_()
{
  // Read in molar masses of NaCl and H2O
  double M_salt;
  if (eos_plist_.isParameter("salt molar mass [kg mol^-1]")) {
    M_salt = eos_plist_.get<double>("salt molar mass [kg mol^-1]");
  } else {
    M_salt = eos_plist_.get<double>("salt molar mass [g mol^-1]", 58.5) * 1e-3;
  }

  double M_water;
  if (eos_plist_.isParameter("water molar mass [kg mol^-1]")) {
    M_water = eos_plist_.get<double>("water molar mass [kg mol^-1]");
  } else {
    M_water = eos_plist_.get<double>("water molar mass [g mol^-1]", 18.0153) * 1e-3;
  }

  // mass density of pure H2O
  rho_f_ = eos_plist_.get<double>("freshwater mass density [kg m^-3]", 1000.0);

  // reference rho and C of saltwater
  double rho_ref = eos_plist_.get<double>("reference saltwater mass density [kg m^-3]", 1025);
  double C_ref = eos_plist_.get<double>("reference saltwater concentration [kg m^-3]", 35.);
  E_ = (rho_ref - rho_f_) / C_ref;
};

} // namespace Relations
} // namespace Amanzi
