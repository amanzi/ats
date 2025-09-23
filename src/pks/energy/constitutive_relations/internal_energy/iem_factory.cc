/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for IEM implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "iem_factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

// method for instantiating IEM implementations
Teuchos::RCP<IEM>
IEMFactory::createIEM(Teuchos::ParameterList& plist)
{
  std::string iem_typename = plist.get<std::string>("IEM type");
  return Teuchos::rcp(CreateInstance(iem_typename, plist));
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
