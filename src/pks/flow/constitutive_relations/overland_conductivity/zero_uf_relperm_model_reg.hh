/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "zero_uf_relperm_model.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<SurfaceRelPermModel, ZeroUFRelPermModel>
  ZeroUFRelPermModel::reg_("unfrozen fraction rel perm, limit zero");

} // namespace Flow
} // namespace Amanzi
