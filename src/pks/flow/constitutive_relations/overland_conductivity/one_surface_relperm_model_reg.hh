/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "one_surface_relperm_model.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<SurfaceRelPermModel, OneSurfaceRelPermModel>
  OneSurfaceRelPermModel::reg_("one surface rel perm");

} // namespace Flow
} // namespace Amanzi
