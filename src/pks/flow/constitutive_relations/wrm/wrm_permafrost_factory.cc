/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

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
  std::string model_typename = plist.get<std::string>("permafrost wrm type");
  Teuchos::RCP<WRMPermafrostModel> model = Teuchos::rcp(CreateInstance(model_typename, plist));
  model->set_WRM(wrm);
  return model;
};

} // namespace Flow
} // namespace Amanzi
