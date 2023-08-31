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

#ifndef AMANZI_UPWINDING_FLUXSPLITDENOMINATOR_SCHEME_
#define AMANZI_UPWINDING_FLUXSPLITDENOMINATOR_SCHEME_

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindFluxSplitDenominator : public Upwinding {
 public:
  UpwindFluxSplitDenominator(const std::string& pkname,
                             const Tag& tag,
                             const std::string& flux,
                             const std::string& slope,
                             const std::string& manning_coef,
                             double flux_epsilon,
                             double slope_regularization);

  virtual void Update(const CompositeVector& cells,
                      CompositeVector& faces,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const override;


  void CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                    const CompositeVector& flux,
                                    const CompositeVector& slope,
                                    const CompositeVector& manning_coef,
                                    CompositeVector& face_coef,
                                    const Teuchos::Ptr<Debugger>& db) const;

  virtual std::string CoefficientLocation() const override { return "upwind: face"; }

 private:
  Tag tag_;
  std::string pkname_;
  std::string flux_;
  std::string slope_;
  std::string manning_coef_;
  double flux_eps_;
  double slope_regularization_;
  std::string ponded_depth_;
};

} // namespace Operators
} // namespace Amanzi

#endif
