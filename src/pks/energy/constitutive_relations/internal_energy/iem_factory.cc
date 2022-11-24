/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for IEM implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "iem_factory.hh"

namespace Amanzi {
namespace Energy {

// method for instantiating IEM implementations
Teuchos::RCP<IEM>
IEMFactory::createIEM(Teuchos::ParameterList& plist)
{
  std::string iem_typename = plist.get<std::string>("IEM type");
  return Teuchos::rcp(CreateInstance(iem_typename, plist));
};

} // namespace Energy
} // namespace Amanzi
