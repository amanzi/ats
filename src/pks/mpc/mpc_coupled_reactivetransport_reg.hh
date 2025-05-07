/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

// ReactiveTransport_PK registration
#include "mpc_coupled_reactivetransport.hh"

namespace Amanzi {

RegisteredPKFactory<MPCCoupledReactiveTransport>
  MPCCoupledReactiveTransport::reg_("surface subsurface reactive transport");

} // namespace Amanzi
