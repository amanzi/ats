/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces. Gives priority to the harmonic average when feasible.
// -----------------------------------------------------------------------------

#ifndef AMANZI_UPWINDING_FLUXHARMONICMEAN_SCHEME_
#define AMANZI_UPWINDING_FLUXHARMONICMEAN_SCHEME_

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindFluxHarmonicMean : public Upwinding {

public:

  UpwindFluxHarmonicMean(std::string pkname,
                         std::string cell_coef,
                         std::string face_coef,
                         std::string flux,
                         double flux_epsilon);
  
  virtual void Update(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<Debugger>& db=Teuchos::null);

  void CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const CompositeVector& flux,
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
  double flux_eps_;
};

} // namespace
} // namespace

#endif
