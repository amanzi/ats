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


template <>
const std::string StrongMPC<PK_BDF_Default>::type = "strong MPC";
template <>
const std::string StrongMPC<PK_PhysicalBDF_Default>::type =
  "physical strong MPC"; // this should not be used?


RegisteredPKFactory<WeakMPC> WeakMPC::reg_("weak MPC");

template <>
RegisteredPKFactory<StrongMPC<PK_BDF_Default>>
  StrongMPC<PK_BDF_Default>::reg_(StrongMPC<PK_BDF_Default>::type);

template <>
RegisteredPKFactory<StrongMPC<PK_PhysicalBDF_Default>>
  StrongMPC<PK_PhysicalBDF_Default>::reg_(StrongMPC<PK_PhysicalBDF_Default>::type);


RegisteredPKFactory<MPCCoupledWater> MPCCoupledWater::reg_("coupled water");

} // namespace Amanzi
