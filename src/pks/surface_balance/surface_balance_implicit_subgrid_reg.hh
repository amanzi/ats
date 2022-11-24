/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Ethan Coon
 *
 * ------------------------------------------------------------------------- */


#include "surface_balance_implicit_subgrid.hh"

namespace Amanzi {
namespace SurfaceBalance {

RegisteredPKFactory<ImplicitSubgrid> ImplicitSubgrid::reg_("surface balance implicit subgrid");

} // namespace SurfaceBalance
} // namespace Amanzi
