/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "mpc_permafrost_split_flux.hh"

namespace Amanzi {

RegisteredPKFactory<MPCPermafrostSplitFlux> MPCPermafrostSplitFlux::reg_(
  "operator split permafrost");

} // namespace Amanzi
