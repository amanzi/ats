/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "wrm_constant.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMConstant> WRMConstant::factory_("constant");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
