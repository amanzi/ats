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

#ifndef AMANZI_UPWINDING_POTENTIALDIFFERENCE_SCHEME_
#define AMANZI_UPWINDING_POTENTIALDIFFERENCE_SCHEME_

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindPotentialDifference : public Upwinding {

public:

  UpwindPotentialDifference(std::string pkname,
                            std::string cell_coef,
                            std::string face_coef,
                            std::string potential,
                            std::string overlap=std::string(""));

  void Update(const Teuchos::Ptr<State>& S,
              const Teuchos::Ptr<Debugger>& db=Teuchos::null);


  void CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const CompositeVector& potential,
        const CompositeVector& overlap,
        CompositeVector& face_coef);

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
  std::string potential_;
  std::string overlap_;
};

} // namespace
} // namespace

#endif
