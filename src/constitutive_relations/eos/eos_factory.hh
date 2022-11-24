/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef PK_FLOW_EOS_FACTORY_HH_
#define PK_FLOW_EOS_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "eos.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Relations {

class EOSFactory : public Amanzi::Utils::Factory<EOS> {
 public:
  Teuchos::RCP<EOS> createEOS(Teuchos::ParameterList& plist);
};

} // namespace Relations
} // namespace Amanzi

#endif
