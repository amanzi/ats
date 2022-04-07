/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

MPCDomainSet registrations.
------------------------------------------------------------------------- */

#include "mpc_domain_set.hh"

namespace Amanzi {

RegisteredPKFactory<MPCDomainSet> MPCDomainSet::reg_("domain set weak MPC");

} // namespace
