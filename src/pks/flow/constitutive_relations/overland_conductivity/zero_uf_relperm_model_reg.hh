/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the conductivity of surface flow.

*/

#include "zero_uf_relperm_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// registry of method
Utils::RegisteredFactory<SurfaceRelPermModel, ZeroUFRelPermModel> ZeroUFRelPermModel::reg_(
  "unfrozen fraction rel perm, limit zero");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
