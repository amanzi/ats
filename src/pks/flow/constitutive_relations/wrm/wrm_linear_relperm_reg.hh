/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  WRM which calls another WRM for saturation but sets 0 rel perm.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_linear_relperm.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMLinearRelPerm> WRMLinearRelPerm::factory_("linear rel perm");

} // namespace Flow
} // namespace Amanzi
