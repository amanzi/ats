/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

// method for instantiating WRM implementations
Teuchos::RCP<WRMPermafrostModel>
WRMPermafrostFactory::createWRMPermafrostModel(Teuchos::ParameterList& plist,
                                               const Teuchos::RCP<WRM>& wrm)
{
  std::string model_typename = plist.get<std::string>("permafrost WRM type");
  Teuchos::RCP<WRMPermafrostModel> model = Teuchos::rcp(CreateInstance(model_typename, plist));
  model->set_WRM(wrm);
  return model;
};

} // namespace Flow
} // namespace Amanzi
