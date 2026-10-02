/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

Accumulates a field pointwise over time, computing either a time-integrated
average, a running maximum, or a running minimum.

This evaluator is designed for use with subcycled MPCs.  It lives at
``tag`` (typically ``Tags::NEXT``) and depends on the source field at
``accumulated_tag`` (the subcycled inner tag, e.g. ``flow_next``).  The
``TimeAdvancer`` driving the inner subcycle loop is responsible for calling
``Reset()`` at the start of each outer timestep and ``Update()`` after each
successful inner step.

For ``integral`` mode the result is the time-weighted mean:

  new_val = (old_val * acc_dt + source * dt) / (acc_dt + dt)

For ``min`` / ``max`` modes the result is the running pointwise
minimum / maximum; accumulated time is still tracked so that ``Reset()``
can initialise correctly.

Both the result ``CompositeVector`` and the ``KEY_accumulated_dt`` scalar are
checkpointed so that the accumulation survives restarts mid outer-step.

`"evaluator type`" = `"time accumulated`"

.. _evaluator-time-accumulated-spec:
.. admonition:: evaluator-time-accumulated-spec

   * `"accumulation type`" ``[string]`` **"integral"** One of ``integral``,
     ``min``, or ``max``.

   * `"accumulated key`" ``[string]`` **my_key** Key of the field to
     accumulate.  Defaults to the same key as this evaluator.

   * `"accumulated tag`" ``[string]`` **required** Tag at which the source
     field is evaluated (the subcycled inner tag).

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

class EvaluatorTimeAccumulated : public EvaluatorSecondary {
 public:
  explicit EvaluatorTimeAccumulated(Teuchos::ParameterList& plist);
  EvaluatorTimeAccumulated(const EvaluatorTimeAccumulated& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  bool IsDifferentiableWRT(const State& S,
                           const Key& wrt_key,
                           const Tag& wrt_tag) const override
  {
    return false;
  }

  void EnsureCompatibility(State& S) override;

  // Zeros the accumulator and clears the lazy-evaluation cache.
  // Called by TimeAdvancer before each outer timestep's inner loop.
  void Reset(State& S);

 protected:
  void Update_(State& S) override;
  void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override {}

 protected:
  Key accumulated_key_;
  Tag accumulated_tag_;
  std::string accumulation_type_;  // "integral", "min", "max"

  // my_keys_[0] = {key, tag}         — the CV result
  // my_keys_[1] = {acc_dt_key, tag}  — the accumulated-dt scalar

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorTimeAccumulated> reg_;
};

} // namespace Relations
} // namespace Amanzi
