/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

// -----------------------------------------------------------------------------
// ATS
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#ifndef AMANZI_UPWINDING_SCHEME_
#define AMANZI_UPWINDING_SCHEME_

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "dbc.hh"
#include "OperatorDefs.hh"
#include "CompositeVector.hh"

namespace Amanzi {

// forward declaration
class State;
class Debugger;
class CompositeVector;

namespace Operators {

enum UpwindMethod {
  UPWIND_METHOD_CENTERED = 0,
  UPWIND_METHOD_GRAVITY,
  UPWIND_METHOD_TOTAL_FLUX,
  UPWIND_METHOD_NO_DENOMINATOR,
  UPWIND_METHOD_ARITHMETIC_MEAN,
  UPWIND_METHOD_POTENTIAL_DIFFERENCE
};

class Upwinding {
 public:
  virtual ~Upwinding() = default;

  virtual void Update(const CompositeVector& data,
                      CompositeVector& uw_data,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const = 0;

  virtual void Update(const CompositeVector& cells,
                      const std::string cell_component,
                      CompositeVector& faces,
                      const std::string face_component,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const
  {
    AMANZI_ASSERT(0);
  }

  virtual void UpdateDerivatives(
    const Teuchos::Ptr<State>& S,
    std::string potential_key,
    const CompositeVector& dconductivity,
    const std::vector<int>& bc_markers,
    const std::vector<double>& bc_values,
    std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>>* Jpp_faces) const
  {
    AMANZI_ASSERT(0);
  }

  virtual std::string CoefficientLocation() const = 0;
};

} // namespace Operators
} // namespace Amanzi

#endif
