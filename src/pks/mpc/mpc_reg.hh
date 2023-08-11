/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/


#include "weak_mpc.hh"
#include "strong_mpc.hh"
#include "mpc_coupled_water.hh"

namespace Amanzi {

RegisteredPKFactory<WeakMPC> WeakMPC::reg_("weak MPC");
template <>
RegisteredPKFactory<StrongMPC<PK_BDF_Default>> StrongMPC<PK_BDF_Default>::reg_("strong MPC");
RegisteredPKFactory<MPCCoupledWater> MPCCoupledWater::reg_("coupled water");

} // namespace Amanzi
