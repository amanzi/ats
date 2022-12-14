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

#include <string>
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {

// method for instantiating WRM implementations
Teuchos::RCP<WRM>
WRMFactory::createWRM(Teuchos::ParameterList& plist)
{
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

} // namespace Flow
} // namespace Amanzi
