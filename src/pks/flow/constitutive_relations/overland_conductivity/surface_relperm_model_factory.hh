/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef PK_FLOW_SURFACE_RELPERM_FACTORY_HH_
#define PK_FLOW_SURFACE_RELPERM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "surface_relperm_model.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class SurfaceRelPermModelFactory : public Utils::Factory<SurfaceRelPermModel> {
 public:
  Teuchos::RCP<SurfaceRelPermModel> createModel(Teuchos::ParameterList& plist);
};

} // namespace Flow
} // namespace Amanzi

#endif
