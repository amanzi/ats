/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Upwinds a field based on the sign of a potential (e.g. pressure) difference.
/*!

Upwind decision is based on the sign of potential_[cells[0]] -
potential_[cells[1]] (the fixed topological cell ordering, not a
flow-dependent orientation), with a smooth transition region near zero
difference sized by an "overlap" field (defaults to the potential itself).

KEYS:

- `"field key`" The field to upwind (on cells).
- `"potential key`" The potential (e.g. pressure) used to decide direction
  (on cells) -- not differentiated through.

*/

#ifndef AMANZI_UPWINDING_EVALUATOR_UPWIND_POTENTIAL_DIFFERENCE_HH_
#define AMANZI_UPWINDING_EVALUATOR_UPWIND_POTENTIAL_DIFFERENCE_HH_

#include "Factory.hh"
#include "evaluator_upwinding.hh"

namespace Amanzi {
namespace Relations {

class EvaluatorUpwindPotentialDifference : public EvaluatorUpwinding {
 public:
  explicit EvaluatorUpwindPotentialDifference(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorUpwindPotentialDifference(const EvaluatorUpwindPotentialDifference& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual std::string getType() const override { return eval_type; }

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override;

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;
  virtual void UpwindDerivative_(const State& S,
                                 const CompositeVector& dfield_dp,
                                 CompositeVector& result) const override;

 protected:
  Key potential_key_;

  static const std::string eval_type;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorUpwindPotentialDifference> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
