/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
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
                            const std::string& face_coef,
                            const std::string& slope,
                            const std::string& manning_coef,
                            const std::string& ponded_depth,
                            const std::string& elevation,
                            const std::string& density,
                            double slope_regularization,
                            double manning_exp);

  virtual void Update(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<Debugger>& db=Teuchos::null);


  void CalculateCoefficientsOnFaces(
        const CompositeVector& slope,
        const CompositeVector& manning_coef,
        const CompositeVector& ponded_depth,
        const CompositeVector& elevation,
        const CompositeVector& density,
        const Teuchos::Ptr<CompositeVector>& face_coef,
        const Teuchos::Ptr<Debugger>& db);

  virtual void
  UpdateDerivatives(const Teuchos::Ptr<State>& S,
                    const std::string& potential_key,
                    const CompositeVector& dconductivity,
                    const std::vector<int>& bc_markers,
                    const std::vector<double>& bc_values,
                    std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > >* Jpp_faces) const;

  virtual std::string
  CoefficientLocation() { return "upwind: face"; }

private:

  std::string pkname_;
  std::string face_coef_;
  std::string slope_;
  std::string manning_coef_;
  std::string ponded_depth_;
  std::string elevation_;
  std::string density_;

  double slope_regularization_;
  double manning_exp_;

};

} // namespace
} // namespace


