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

#ifndef PK_FLOW_WRM_FACTORY_HH_
#define PK_FLOW_WRM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class WRMFactory : public Utils::Factory<WRM> {
 public:
  Teuchos::RCP<WRM> createWRM(Teuchos::ParameterList& plist);
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
