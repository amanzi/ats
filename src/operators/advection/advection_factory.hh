/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

#ifndef OPERATOR_ADVECTION_ADVECTION_FACTORY_HH_
#define OPERATOR_ADVECTION_ADVECTION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "CompositeVector.hh"

#include "advection.hh"

namespace Amanzi {
namespace Operators {

class AdvectionFactory {
 public:
  Teuchos::RCP<Advection>
  create(Teuchos::ParameterList& advect_plist, const Teuchos::RCP<const AmanziMesh::Mesh> mesh);
};

} // namespace Operators
} // namespace Amanzi

#endif
