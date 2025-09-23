/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  Constant density/viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

*/

#include "eos_constant.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

EOSConstant::EOSConstant(Teuchos::ParameterList& eos_plist)
  : eos_plist_(eos_plist)
{
  InitializeFromPlist_();
};

void
EOSConstant::InitializeFromPlist_()
{
  // defaults to water
  if (eos_plist_.isParameter("molar mass [kg mol^-1]")) {
    M_ = eos_plist_.get<double>("molar mass [kg mol^-1]");
  } else {
    M_ = eos_plist_.get<double>("molar mass [g mol^-1]", 18.0153) * 1.e-3;
  }

  if (eos_plist_.isParameter("density [mol m^-3]")) {
    rho_ = eos_plist_.get<double>("density [mol m^-3]") * M_;
  } else {
    rho_ = eos_plist_.get<double>("density [kg m^-3]");
  }
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
