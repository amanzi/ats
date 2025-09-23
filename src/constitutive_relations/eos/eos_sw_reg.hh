/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  ATS

  EOS for an salt water (does not implement viscosity at this point!)

*/

#include "eos_factory.hh"
#include "eos_sw.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

Utils::RegisteredFactory<EOS, EOS_SW> EOS_SW::factory_("salt water");

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
