/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#ifndef PK_FLOW_WRM_FACTORY_HH_
#define PK_FLOW_WRM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMFactory : public Utils::Factory<WRM> {
 public:
  Teuchos::RCP<WRM> createWRM(Teuchos::ParameterList& plist);
};

} // namespace Flow
} // namespace Amanzi

#endif
