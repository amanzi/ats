/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for overland flow.
------------------------------------------------------------------------- */

#include "energy_surface_ice.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<EnergySurfaceIce> EnergySurfaceIce::reg_("surface energy");

} // namespace Energy
} // namespace Amanzi
