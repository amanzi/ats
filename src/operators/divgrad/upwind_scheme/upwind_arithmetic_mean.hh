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

#ifndef AMANZI_UPWINDING_ARITHMETICMEAN_SCHEME_
#define AMANZI_UPWINDING_ARITHMETICMEAN_SCHEME_

#include "Key.hh"
#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindArithmeticMean : public Upwinding {

public:

  UpwindArithmeticMean(Key pkname,
                     Key cell_coef,
                     Key face_coef);
  UpwindArithmeticMean(const UpwindArithmeticMean& other) = delete;
  UpwindArithmeticMean& operator=(const UpwindArithmeticMean& other) = delete;

  virtual void Update(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<Debugger>& db=Teuchos::null);

  void CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        CompositeVector& face_coef);

  virtual void
  UpdateDerivatives(const Teuchos::Ptr<State>& S, 
                    Key potential_key,
                    const CompositeVector& dconductivity,
                    const std::vector<int>& bc_markers,
                    const std::vector<double>& bc_values,
                    std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > >* Jpp_faces) const;

  virtual std::string
  CoefficientLocation() { return "upwind: face"; }
  

private:
  
  Key pkname_;
  Key cell_coef_;
  Key face_coef_;
};

} // namespace
} // namespace

#endif
