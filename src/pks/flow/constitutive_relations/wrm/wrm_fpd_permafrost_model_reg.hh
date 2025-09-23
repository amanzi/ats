/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Painter's permafrost model with freezing point depression.

 */

#include "wrm.hh"
#include "wrm_fpd_permafrost_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {


// registry of method
Utils::RegisteredFactory<WRMPermafrostModel, WRMFPDPermafrostModel> WRMFPDPermafrostModel::factory_(
  "fpd permafrost model");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
