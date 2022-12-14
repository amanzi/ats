/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Constant density/viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_constant.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<EOS, EOSConstant> EOSConstant::factory_("constant");

} // namespace Relations
} // namespace Amanzi
