/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Linear density/viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_linear.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<EOS, EOSLinear> EOSLinear::factory_("linear");

} // namespace Relations
} // namespace Amanzi
