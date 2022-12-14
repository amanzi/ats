/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "one_uf_relperm_model.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<SurfaceRelPermModel, OneUFRelPermModel>
  OneUFRelPermModel::reg_("unfrozen fraction rel perm, limit one");

} // namespace Flow
} // namespace Amanzi
