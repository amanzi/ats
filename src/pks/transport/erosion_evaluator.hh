/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  The erosion evaluator gets the erosion rates.


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
  explicit ErosionRateEvaluator(Teuchos::ParameterList& plist);
  ErosionRateEvaluator(const ErosionRateEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const override;

  // virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
  //         const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;
  // virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);
  // virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S){};

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  double tau_e_;
  double Qe_0_;
  double gamma_;
  double lambda_, umax_, xi_;
  double Cf_;

  Key velocity_key_;

  static Utils::RegisteredFactory<Evaluator, ErosionRateEvaluator> factory_;
};

} // namespace Amanzi

#endif
