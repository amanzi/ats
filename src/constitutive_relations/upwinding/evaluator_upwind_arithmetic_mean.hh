/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Upwinds a field as the simple arithmetic mean of its two neighbor cells.
/*!

No directionality -- the face value is always 0.5*(cell0 + cell1), so both
derivative slots are always nonzero (0.5 each), unlike the hard-upwind rules.

KEYS:

- `"field key`" The field to upwind (on cells).

*/

#ifndef AMANZI_UPWINDING_EVALUATOR_UPWIND_ARITHMETIC_MEAN_HH_
#define AMANZI_UPWINDING_EVALUATOR_UPWIND_ARITHMETIC_MEAN_HH_

#include "Factory.hh"
#include "evaluator_upwinding.hh"

namespace Amanzi {
namespace Relations {

class EvaluatorUpwindArithmeticMean : public EvaluatorUpwinding {
 public:
  explicit EvaluatorUpwindArithmeticMean(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorUpwindArithmeticMean(const EvaluatorUpwindArithmeticMean& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual std::string getType() const override { return eval_type; }

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;
  virtual void UpwindDerivative_(const State& S,
                                 const CompositeVector& dfield_dp,
                                 CompositeVector& result) const override;

 protected:
  static const std::string eval_type;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorUpwindArithmeticMean> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
