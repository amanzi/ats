/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGYRELATIONS_IEM_FACTORY_
#define AMANZI_ENERGYRELATIONS_IEM_FACTORY_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "iem.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class IEMFactory : public Utils::Factory<IEM> {
 public:
  Teuchos::RCP<IEM> createIEM(Teuchos::ParameterList& plist);
};

} // namespace Energy
} // namespace Amanzi

#endif
