/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

Computes a maximum (or minimum) of a field, over time in a pointwise way.

.. _evaluator-time-max-spec:
.. admonition:: evaluator-time-max-spec

   * `"operator`" ``[string]`` **max** One of max or min.

   KEYS:

   - `"dependency key`" **default** Defaults to the same domain/variable as
     this key, without a "max_" or "min_" prefix
*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class TimeMaxEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit TimeMaxEvaluator(Teuchos::ParameterList& plist);
  TimeMaxEvaluator(const TimeMaxEvaluator& other) = default;
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

 protected:
  std::string operator_;
  bool evaluated_once_;

 private:
  static Utils::RegisteredFactory<Evaluator, TimeMaxEvaluator> reg_;
};

} // namespace Relations
} // namespace Amanzi
