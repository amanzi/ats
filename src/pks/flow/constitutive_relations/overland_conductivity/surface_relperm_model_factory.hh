/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef PK_FLOW_SURFACE_RELPERM_FACTORY_HH_
#define PK_FLOW_SURFACE_RELPERM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "surface_relperm_model.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class SurfaceRelPermModelFactory : public Utils::Factory<SurfaceRelPermModel> {
 public:
  Teuchos::RCP<SurfaceRelPermModel> createModel(Teuchos::ParameterList& plist);
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
