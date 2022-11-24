/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Ethan Coon, Adam Atchley, Satish Karra
 *
 * ------------------------------------------------------------------------- */


#include "surface_balance_base.hh"

namespace Amanzi {
namespace SurfaceBalance {

RegisteredPKFactory<SurfaceBalanceBase> SurfaceBalanceBase::reg_("general surface balance");

} // namespace SurfaceBalance
} // namespace Amanzi
