/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#ifndef AMANZI_UPWINDING_CELLCENTERED_SCHEME_
#define AMANZI_UPWINDING_CELLCENTERED_SCHEME_

#include "Key.hh"
#include "Tag.hh"
#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindCellCentered : public Upwinding {
 public:
  UpwindCellCentered(const std::string& pkname, const Tag& tag);

  virtual void Update(const CompositeVector& cells,
                      CompositeVector& faces,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const override;

  virtual std::string CoefficientLocation() const override { return "standard: cell"; }

 private:
  std::string pkname_;
  Tag tag_;
};

} // namespace Operators
} // namespace Amanzi

#endif
