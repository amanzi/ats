/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: ...
           Ethan Coon (ATS version) (coonet@ornl.gov)
*/

#ifndef AMANZI_TEST_PK_BC_FACTORY_HH_
#define AMANZI_TEST_PK_BC_FACTORY_HH_


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "bc_factory.hh"

namespace Amanzi {
namespace TestPKs {

class TestPKBCFactory : public BCFactory {
 public:
  TestPKBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                  const Teuchos::ParameterList& plist)
    : BCFactory(mesh, plist)
  {}

  Teuchos::RCP<Functions::BoundaryFunction> CreateDirichlet() const
  {
    return CreateWithFunction("dirichlet", "boundary data");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateNeumann() const
  {
    return CreateWithFunction("neumann", "outward flux");
  }
};

} // namespace TestPKs
} // namespace Amanzi

#endif // AMANZI_FLOW_BC_FACTORY_HH_
