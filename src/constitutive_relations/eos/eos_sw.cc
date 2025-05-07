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
  double mol_ratio = params[0];
  return rho_f_ + E_ * mol_ratio;
};

double
EOS_SW::DMassDensityDMolarRatio(std::vector<double>& params)
{
  return E_;
};

double
EOS_SW::MolarDensity(std::vector<double>& params) // H2O
{
  //return MassDensity(params) / (M_water_ * (1 - C) + M_salt_ * C);
  return n_l_;
};

double
EOS_SW::DMolarDensityDMolarRatio(std::vector<double>& params)
{
  return 0.;
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

  // mass density of H2O
  rho_f_ = eos_plist_.get<double>("fresh water mass density [kg m^-3]", 1000.0);

  // molar density of H2O
  n_l_ = rho_f_ / M_water;

  // reference rho and C of saltwater
  double rho_max = eos_plist_.get<double>("reference saltwater mass density [kg m^-3]", 1025);
  double C_max = eos_plist_.get<double>("reference salinity [kg salt m^-3]", 35.);

  // convert to mole ratio
  double C_max_mol_m3 = C_max / M_salt; // [mol NaCl m^-3]
  double mol_ratio_max = C_max_mol_m3 * n_l_; // [mol NaCl (mol H2O)^-1]
  E_ = (rho_max - rho_f_) / mol_ratio_max; // [kg m^-3 mol H2O (mol NaCl)^-1]
};

} // namespace Relations
} // namespace Amanzi
