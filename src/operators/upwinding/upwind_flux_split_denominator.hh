/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
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

  UpwindFluxSplitDenominator(std::string pkname,
                             std::string cell_coef,
                             std::string face_coef,
                             std::string flux,
                             double flux_epsilon,
                             std::string slope,
                             std::string manning_coef,
                             double slope_regularization);

  virtual void Update(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<Debugger>& db=Teuchos::null);


  void CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const CompositeVector& flux,
        const CompositeVector& slope,
        const CompositeVector& manning_coef,
        const Teuchos::Ptr<CompositeVector>& face_coef,
        const Teuchos::Ptr<Debugger>& db);

  virtual void
  UpdateDerivatives(const Teuchos::Ptr<State>& S, 
                    std::string potential_key,
                    const CompositeVector& dconductivity,
                    const std::vector<int>& bc_markers,
                    const std::vector<double>& bc_values,
                    std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > >* Jpp_faces) const;

  virtual std::string
  CoefficientLocation() { return "upwind: face"; }
  
private:

  std::string pkname_;
  std::string cell_coef_;
  std::string face_coef_;
  std::string flux_;
  double flux_eps_;
  std::string slope_;
  std::string manning_coef_;
  double slope_regularization_;
  std::string ponded_depth_;
};

} // namespace
} // namespace

#endif
