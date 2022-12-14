
#include "wrm.hh"
#include "wrm_implicit_permafrost_model.hh"


namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<WRMPermafrostModel, WRMImplicitPermafrostModel>
  WRMImplicitPermafrostModel::factory_("permafrost model");

} // namespace Flow
} // namespace Amanzi
