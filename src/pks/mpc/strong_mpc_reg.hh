/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived StrongMPC class.  Is both a PK and a Model
Evalulator, providing needed methods for BDF time integration of the coupled
system.

Completely automated and generic to any sub PKs, this uses a block diagonal
preconditioner.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */


#include "strong_mpc.hh"

namespace Amanzi {

template <>
RegisteredPKFactory<StrongMPC<PK_BDF_Default>> StrongMPC<PK_BDF_Default>::reg_("strong MPC");

} // namespace Amanzi
