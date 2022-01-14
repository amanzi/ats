/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.


  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_EROSIONRATE_EVALUATOR_
#define AMANZI_EROSIONRATE_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class ErosionRateEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  ErosionRateEvaluator(Teuchos::ParameterList& plist);

  ErosionRateEvaluator(const ErosionRateEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;
  
  // virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
  //         const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;

  // virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);

  //virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S){};

protected:

    // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key,
                                               const Teuchos::Ptr<CompositeVector>& result);

  double tau_e_;
  double Qe_0_;
  double gamma_;
  double lambda_, umax_, xi_;

  Key velocity_key_;

  static Utils::RegisteredFactory<Evaluator,ErosionRateEvaluator> factory_;

};

} //namespace

#endif
