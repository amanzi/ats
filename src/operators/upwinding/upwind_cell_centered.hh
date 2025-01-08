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
