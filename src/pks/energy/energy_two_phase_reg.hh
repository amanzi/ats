/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#include "energy_two_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<TwoPhase> TwoPhase::reg_("two-phase energy");

} // namespace Energy
} // namespace Amanzi
