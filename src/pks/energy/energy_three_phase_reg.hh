/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "energy_three_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<ThreePhase> ThreePhase::reg_("three-phase energy");

} // namespace Energy
} // namespace Amanzi
