/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "energy_interfrost.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<InterfrostEnergy> InterfrostEnergy::reg_("interfrost energy");

} // namespace Energy
} // namespace Amanzi
