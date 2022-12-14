/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "mpc_subcycled.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSubcycled> MPCSubcycled::reg_("subcycling MPC");

} // namespace Amanzi
