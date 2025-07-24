/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*

McKenzie et al. (2007)'s soil freezing curve

 */

#include "wrm.hh"
#include "wrm_mck_permafrost_model.hh"

namespace Amanzi {
namespace Flow {


// registry of method
Utils::RegisteredFactory<WRMPermafrostModel, WRMMCKPermafrostModel> WRMMCKPermafrostModel::factory_(
  "mck permafrost model");

} // namespace Flow
} // namespace Amanzi
