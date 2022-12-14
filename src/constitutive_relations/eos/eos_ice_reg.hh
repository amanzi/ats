/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid ice.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_ice.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<EOS, EOSIce> EOSIce::factory_("ice");

} // namespace Relations
} // namespace Amanzi
