/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  Constant viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

*/

#include "viscosity_constant.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<ViscosityRelation, ViscosityConstant> ViscosityConstant::factory_(
  "constant");

} // namespace Relations
} // namespace Amanzi
