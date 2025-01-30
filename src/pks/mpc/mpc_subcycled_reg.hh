/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Implementation for the derived WeakMPC class.  Provides only the advance()
method missing from MPC.hh.  In weak coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "mpc_subcycled.hh"
#include "mpc_flow_transport.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSubcycled> MPCSubcycled::reg_("subcycling MPC");
RegisteredPKFactory<MPCFlowTransport> MPCFlowTransport::reg_("coupled flow and transport");

} // namespace Amanzi
