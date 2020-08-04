/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

DomainSetMPC registrations.
------------------------------------------------------------------------- */

#include "DomainSetMPC.hh"

namespace Amanzi {

RegisteredPKFactory<DomainSetMPC> DomainSetMPC::reg_("domain set weak MPC");

} // namespace
