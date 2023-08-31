/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  Linear density/viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

*/

#include "eos_linear.hh"

namespace Amanzi {
namespace Relations {

EOSLinear::EOSLinear(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist)
{
  InitializeFromPlist_();
};

void
EOSLinear::InitializeFromPlist_()
{
  // defaults to water
  if (eos_plist_.isParameter("molar mass [kg/mol]")) {
    M_ = eos_plist_.get<double>("molar mass [kg/mol]");
  } else {
    M_ = eos_plist_.get<double>("molar mass [g/mol]", 18.0153) * 1.e-3;
  }

  if (eos_plist_.isParameter("density [mol/m^3]")) {
    rho_ = eos_plist_.get<double>("density [mol/m^3]") * M_;
  } else {
    rho_ = eos_plist_.get<double>("density [kg/m^3]");
  }

  beta_ = eos_plist_.get<double>("compressibility [1/Pa]");
};

} // namespace Relations
} // namespace Amanzi
