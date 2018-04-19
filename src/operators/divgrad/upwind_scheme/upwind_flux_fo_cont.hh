/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
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

  UpwindFluxFOCont(std::string pkname,
                             std::string cell_coef,
                             std::string face_coef,
                             std::string flux,
                             std::string slope,
                             std::string manning_coef,
                             std::string elevation,
                             double slope_regularization,
                             double manning_exp);

  virtual void Update(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<Debugger>& db=Teuchos::null);


  void CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const CompositeVector& flux,
        const CompositeVector& slope,
        const CompositeVector& manning_coef,
        const CompositeVector& elevation,
        CompositeVector& face_coef,
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
  std::string slope_;
  std::string manning_coef_;
  std::string elevation_;
  double slope_regularization_;
  double manning_exp_;
};

} // namespace
} // namespace

#endif
