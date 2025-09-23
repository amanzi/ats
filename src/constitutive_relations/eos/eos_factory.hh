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

#ifndef PK_FLOW_EOS_FACTORY_HH_
#define PK_FLOW_EOS_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "eos.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

class EOSFactory : public Amanzi::Utils::Factory<EOS> {
 public:
  Teuchos::RCP<EOS> createEOS(Teuchos::ParameterList& plist);
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
