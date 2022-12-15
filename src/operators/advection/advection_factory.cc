/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Interface for a general-purpose advection factory.
   ------------------------------------------------------------------------- */

#include "errors.hh"
#include "advection_donor_upwind.hh"
#include "advection_factory.hh"

namespace Amanzi {
namespace Operators {

Teuchos::RCP<Advection>
AdvectionFactory::create(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  std::string method = plist.get<std::string>("Advection method");

  if (method == "donor upwind") {
    return Teuchos::rcp(new AdvectionDonorUpwind(plist, mesh));
  } else {
    std::stringstream emsg;
    emsg << "AdvectionFactory: unknown advection type " << method;
    Errors::Message message(emsg.str());
    Exceptions::amanzi_throw(message);
  }

  return Teuchos::null;
};

} // namespace Operators
} // namespace Amanzi
