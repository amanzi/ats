/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "icy_overland.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

RegisteredPKFactory<IcyOverlandFlow> IcyOverlandFlow::reg_("overland flow with ice");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
