/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for WRM_PERMAFROST implementations.
   ------------------------------------------------------------------------- */

#ifndef _PK_FLOW_WRM_PERMAFROST_FACTORY_HH_
#define _PK_FLOW_WRM_PERMAFROST_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "wrm_permafrost_model.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMPermafrostFactory : public Utils::Factory<WRMPermafrostModel> {
 public:
  Teuchos::RCP<WRMPermafrostModel> createWRMPermafrostModel(Teuchos::ParameterList& plist,
                                                            const Teuchos::RCP<WRM>& wrm);
};

} // namespace Flow
} // namespace Amanzi

#endif
