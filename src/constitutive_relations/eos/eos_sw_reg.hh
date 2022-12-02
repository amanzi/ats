/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

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
