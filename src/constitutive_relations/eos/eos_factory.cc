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

#include <string>
#include "eos_factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

// method for instantiating EOS implementations
Teuchos::RCP<EOS>
EOSFactory::createEOS(Teuchos::ParameterList& plist)
{
  std::string eos_typename = plist.get<std::string>("EOS type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};


} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
