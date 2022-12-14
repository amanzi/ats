/*
Author: Ethan Coon

Painter's permafrost model with freezing point depression.

 */

#include "wrm.hh"
#include "wrm_interfrost_permafrost_model.hh"

namespace Amanzi {
namespace Flow {


// registry of method
Utils::RegisteredFactory<WRMPermafrostModel, WRMInterfrostPermafrostModel>
  WRMInterfrostPermafrostModel::factory_("interfrost permafrost model");

} // namespace Flow
} // namespace Amanzi
