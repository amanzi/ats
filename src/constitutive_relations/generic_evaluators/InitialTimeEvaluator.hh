/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

Stores the value of a field at the initial time.  Note this must be
checkpointed because it can never be reconstructed.

.. _initial_time_evaluator-spec:
.. admonition:: initial_time_evaluator-spec

   * `"initial time`" ``[double]`` **0** Time value to save, in units as below.
   * `"initial time units`" ``[string]`` **s** Units of initial time above.

   KEYS:

   - `"dependency key`" **default** Defaults to the same domain/variable as
     this key, without an "init_" prefix
*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class InitialTimeEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit InitialTimeEvaluator(Teuchos::ParameterList& plist);
  InitialTimeEvaluator(const InitialTimeEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override
  {}

  // must be checkpointed
  virtual void EnsureCompatibility_Flags_(State& S) override;

 protected:
  std::string operator_;
  bool evaluated_once_;
  double time_;

 private:
  static Utils::RegisteredFactory<Evaluator, InitialTimeEvaluator> reg_;
};

} // namespace Relations
} // namespace Amanzi
