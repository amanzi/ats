/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon, Adam Atchley, Satish Karra
*/

#include "surface_balance_base.hh"

namespace Amanzi {
namespace SurfaceBalance {

RegisteredPKFactory<SurfaceBalanceBase> SurfaceBalanceBase::reg_("general surface balance");

} // namespace SurfaceBalance
} // namespace Amanzi
