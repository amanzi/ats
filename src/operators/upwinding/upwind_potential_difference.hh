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

#ifndef AMANZI_UPWINDING_POTENTIALDIFFERENCE_SCHEME_
#define AMANZI_UPWINDING_POTENTIALDIFFERENCE_SCHEME_

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindPotentialDifference : public Upwinding {
 public:
  UpwindPotentialDifference(const std::string& pkname,
                            const Tag& tag,
                            const std::string& potential,
                            const std::string& overlap = std::string(""));

  virtual void Update(const CompositeVector& data,
                      CompositeVector& uw_data,
                      const State& S,
                      const Teuchos::Ptr<Debugger>& db = Teuchos::null) const override;

  void CalculateCoefficientsOnFaces(const CompositeVector& cell_coef,
                                    const CompositeVector& potential,
                                    const CompositeVector& overlap,
                                    CompositeVector& face_coef) const;

  virtual void UpdateDerivatives(
    const Teuchos::Ptr<State>& S,
    std::string potential_key,
    const CompositeVector& dconductivity,
    const std::vector<int>& bc_markers,
    const std::vector<double>& bc_values,
    std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>>* Jpp_faces) const override;

  virtual std::string CoefficientLocation() const override { return "upwind: face"; }

 private:
  std::string pkname_;
  Tag tag_;
  std::string potential_;
  std::string overlap_;
};

} // namespace Operators
} // namespace Amanzi

#endif
