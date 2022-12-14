/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an salt water (does not implement viscosity at this point!)

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_sw.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<EOS, EOS_SW> EOS_SW::factory_("salt water");

} // namespace Relations
} // namespace Amanzi
