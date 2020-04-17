/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Ethan Coon, Adam Atchley, Satish Karra
 *
 * ------------------------------------------------------------------------- */


#include "surface_balance_base.hh"

template<>
Amanzi::RegisteredPKFactory<ATS::SurfaceBalance::SurfaceBalanceBase_Implicit> ATS::SurfaceBalance::SurfaceBalanceBase_Implicit::reg_("general surface balance, implicit");
template<>
Amanzi::RegisteredPKFactory<ATS::SurfaceBalance::SurfaceBalanceBase_Explicit> ATS::SurfaceBalance::SurfaceBalanceBase_Explicit::reg_("general surface balance, explicit");
