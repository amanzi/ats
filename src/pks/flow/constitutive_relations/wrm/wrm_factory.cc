/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {

// method for instantiating WRM implementations
Teuchos::RCP<WRM> WRMFactory::createWRM(Teuchos::ParameterList& plist) {
  std::string wrm_typename;
  // need to deprecate "WRM Type"
  if (plist.isParameter("wrm type"))
    wrm_typename = plist.get<std::string>("wrm type");
  else if (plist.isParameter("WRM Type"))
    wrm_typename = plist.get<std::string>("WRM Type");
  else if (plist.isParameter("WRM type"))
    wrm_typename = plist.get<std::string>("WRM type");
  else
    wrm_typename = plist.get<std::string>("wrm type");
  return Teuchos::rcp(CreateInstance(wrm_typename, plist));
};

} // namespace
} // namespace

