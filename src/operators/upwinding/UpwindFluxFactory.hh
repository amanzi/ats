// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
//
// Factory for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#ifndef AMANZI_UPWINDING_FACTORY_HH_
#define AMANZI_UPWINDING_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "upwinding.hh"

namespace Amanzi {
namespace Operators {

class UpwindFluxFactory {
 public:
  UpwindFluxFactory() {};
  ~UpwindFluxFactory() {};

  Teuchos::RCP<Upwinding> Create(Teuchos::ParameterList& oplist,
          const Teuchos::Ptr<State>& S,
          std::string pkname,
          std::string cell_coef,
          std::string face_coef,
          std::string flux);

};

}  // namespace Operators
}  // namespace Amanzi

#endif
