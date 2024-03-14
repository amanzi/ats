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

#ifndef AMANZI_UPWINDING_FLUXFOCONT_SCHEME_
#define AMANZI_UPWINDING_FLUXFOCONT_SCHEME_

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindFluxFOCont : public Upwinding {
 public:
  UpwindFluxFOCont(const std::string& pkname,
                   const Tag& tag,
                   const Key& flux,
                   const Key& slope,
                   const Key& manning_coef,
                   const Key& elevation,
                   double slope_regularization,
                   double manning_exp);

  virtual void Update(const CompositeVector& data,
                      CompositeVector& uw_data,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const override;

  void CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                    const CompositeVector& flux,
                                    const CompositeVector& slope,
                                    const CompositeVector& manning_coef,
                                    const CompositeVector& elevation,
                                    CompositeVector& face_coef,
                                    const Teuchos::Ptr<Debugger>& db) const;

  virtual std::string CoefficientLocation() const override { return "upwind: face"; }

 private:
  std::string pkname_;
  Tag tag_;
  std::string flux_;
  std::string slope_;
  std::string manning_coef_;
  std::string elevation_;
  double slope_regularization_;
  double manning_exp_;
};

} // namespace Operators
} // namespace Amanzi

#endif
