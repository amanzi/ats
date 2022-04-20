/*

McKenzie et al. (2007)'s soil freezing curve 

 */

#include "wrm.hh"
#include "wrm_mck_permafrost_model.hh"

namespace Amanzi {
namespace Flow {


// registry of method
Utils::RegisteredFactory<WRMPermafrostModel,WRMMCKPermafrostModel> WRMMCKPermafrostModel::factory_("mck permafrost model");

} // namespace Flow
} // namespace Flow
