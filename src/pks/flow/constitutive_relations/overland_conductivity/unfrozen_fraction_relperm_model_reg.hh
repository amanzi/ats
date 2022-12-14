/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_fraction_relperm_model.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<SurfaceRelPermModel, UnfrozenFractionRelPermModel>
  UnfrozenFractionRelPermModel::reg_("unfrozen fraction rel perm");

} // namespace Flow
} // namespace Amanzi
