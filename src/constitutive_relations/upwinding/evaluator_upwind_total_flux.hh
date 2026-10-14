/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Upwinds a field based on the sign of a precomputed flux-direction field.
/*!

The most commonly used upwinding rule.  Face coefficient is set to whichever
neighboring cell's value is upwind of the flux, with a smooth (linear) blend
within a tolerance band around zero flux.

KEYS:

- `"field key`" The field to upwind (on cells).
- `"flux direction key`" The (already gravity-inclusive, if relevant)
  flux-direction field on faces, used only to decide direction -- not
  differentiated through.

*/

#ifndef AMANZI_UPWINDING_EVALUATOR_UPWIND_TOTAL_FLUX_HH_
#define AMANZI_UPWINDING_EVALUATOR_UPWIND_TOTAL_FLUX_HH_

#include "Factory.hh"
#include "evaluator_upwinding.hh"

namespace Amanzi {
namespace Relations {

class EvaluatorUpwindTotalFlux : public EvaluatorUpwinding {
 public:
  explicit EvaluatorUpwindTotalFlux(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorUpwindTotalFlux(const EvaluatorUpwindTotalFlux& other) = default;
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
  Key flux_dir_key_;
  double flux_eps_;

  static const std::string eval_type;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorUpwindTotalFlux> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
