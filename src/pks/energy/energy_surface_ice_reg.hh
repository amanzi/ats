/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Process kernel for energy equation for overland flow.
------------------------------------------------------------------------- */

#include "energy_surface_ice.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

RegisteredPKFactory<EnergySurfaceIce> EnergySurfaceIce::reg_("surface energy");

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
