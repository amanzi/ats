/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

// -----------------------------------------------------------------------------
// ATS
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#pragma once

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindElevationStabilized : public Upwinding {
 public:
  UpwindElevationStabilized(const std::string& pkname,
                            const Tag& tag,
                            const std::string& slope,
                            const std::string& manning_coef,
                            const std::string& ponded_depth,
                            const std::string& elevation,
                            const std::string& density,
                            double slope_regularization,
                            double manning_exp);

  virtual void Update(const CompositeVector& cells,
                      CompositeVector& faces,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const override;

  void CalculateCoefficientsOnFaces(const CompositeVector& slope,
                                    const CompositeVector& manning_coef,
                                    const CompositeVector& ponded_depth,
                                    const CompositeVector& elevation,
                                    const CompositeVector& density,
                                    CompositeVector& face_coef,
                                    const Teuchos::Ptr<Debugger>& db) const;

  virtual std::string CoefficientLocation() const override { return "upwind: face"; }

 private:
  Tag tag_;
  std::string pkname_;
  std::string slope_;
  std::string manning_coef_;
  std::string ponded_depth_;
  std::string elevation_;
  std::string density_;

  double slope_regularization_;
  double manning_exp_;
};

} // namespace Operators
} // namespace Amanzi
