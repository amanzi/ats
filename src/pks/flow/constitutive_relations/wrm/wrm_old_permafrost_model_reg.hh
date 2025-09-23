/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Painter's permafrost model.

 */

#include "wrm.hh"
#include "wrm_old_permafrost_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

Utils::RegisteredFactory<WRMPermafrostModel, WRMOldPermafrostModel> WRMOldPermafrostModel::factory_(
  "old permafrost model");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
