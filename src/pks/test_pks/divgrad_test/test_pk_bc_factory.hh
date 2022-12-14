#ifndef AMANZI_TEST_PK_BC_FACTORY_HH_
#define AMANZI_TEST_PK_BC_FACTORY_HH_

/* -------------------------------------------------------------------------
ATS

Author: ...
    Ethan Coon (ATS version) (ecoon@lanl.gov)

*/

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
