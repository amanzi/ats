/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Constant viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "viscosity_constant.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<ViscosityRelation, ViscosityConstant>
  ViscosityConstant::factory_("constant");

} // namespace Relations
} // namespace Amanzi
