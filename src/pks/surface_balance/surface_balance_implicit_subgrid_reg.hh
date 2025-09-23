/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "surface_balance_implicit_subgrid.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace SurfaceBalance {

RegisteredPKFactory<ImplicitSubgrid> ImplicitSubgrid::reg_("surface balance implicit subgrid");

} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
