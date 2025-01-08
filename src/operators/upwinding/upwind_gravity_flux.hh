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

#ifndef AMANZI_UPWINDING_GRAVITYFLUX_SCHEME_
#define AMANZI_UPWINDING_GRAVITYFLUX_SCHEME_

#include "Epetra_Vector.h"
#include "Tensor.hh"

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindGravityFlux : public Upwinding {
 public:
  UpwindGravityFlux(const std::string& pkname,
                    const Tag& tag,
                    const Teuchos::RCP<std::vector<WhetStone::Tensor>> K);

  virtual void Update(const CompositeVector& cells,
                      CompositeVector& faces,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const override;

  virtual void Update(const CompositeVector& cells,
                      const std::string cell_component,
                      CompositeVector& faces,
                      const std::string face_component,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const override;

  void CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                    const std::string cell_component,
                                    const AmanziGeometry::Point& gravity,
                                    CompositeVector& face_coef,
                                    const std::string face_component) const;

  virtual std::string CoefficientLocation() const override { return "upwind: face"; }

 private:
  std::string pkname_;
  Tag tag_;
  Teuchos::RCP<std::vector<WhetStone::Tensor>> K_;
};

} // namespace Operators
} // namespace Amanzi

#endif
