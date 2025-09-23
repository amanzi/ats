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

#ifndef AMANZI_ENERGYRELATIONS_IEM_FACTORY_
#define AMANZI_ENERGYRELATIONS_IEM_FACTORY_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "iem.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class IEMFactory : public Utils::Factory<IEM> {
 public:
  Teuchos::RCP<IEM> createIEM(Teuchos::ParameterList& plist);
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
