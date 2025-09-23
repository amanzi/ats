/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_vapor_in_gas.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

EOSVaporInGas::EOSVaporInGas(Teuchos::ParameterList& eos_plist)
  : eos_plist_(eos_plist)
{
  InitializeFromPlist_();
}

double
EOSVaporInGas::MolarDensity(std::vector<double>& params)
{
  return gas_eos_->MolarDensity(params);
};

double
EOSVaporInGas::DMolarDensityDT(std::vector<double>& params)
{
  return gas_eos_->DMolarDensityDT(params);
};

double
EOSVaporInGas::DMolarDensityDp(std::vector<double>& params)
{
  return gas_eos_->DMolarDensityDp(params);
};


void
EOSVaporInGas::InitializeFromPlist_()
{
  Teuchos::ParameterList gas_eos_plist = eos_plist_.sublist("gas EOS parameters");
  EOSFactory eos_factory;
  gas_eos_ = eos_factory.createEOS(gas_eos_plist);
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
