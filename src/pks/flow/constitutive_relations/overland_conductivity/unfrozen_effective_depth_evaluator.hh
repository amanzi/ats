/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates the unfrozen mobile depth.
/*!

In freezing conditions, water is only mobile if it is unfrozen.  This evaluator
determines how much water is allowed to flow given that it is partially frozen.

.. math:

   \delta_{mobile} = \delta \chi^{\alpha}

Given a ponded depth, an unfrozen fraction, and an optional power-law exponent,
which we call the ice retardation exponent.

.. _unfrozen-effective-depth-evaluator-spec:
.. admonition:: unfrozen-effective-depth-evaluator-spec

  * `"ice retardation exponent [-]`" ``[double]`` **1.0** exponent alpha
    controlling how quickly ice turns off flow.

  DEPENDENCIES:
  - `"depth`" **DOMAIN-depth**
  - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**

*/


#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class UnfrozenEffectiveDepthModel;

class UnfrozenEffectiveDepthEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit UnfrozenEffectiveDepthEvaluator(Teuchos::ParameterList& plist);
  UnfrozenEffectiveDepthEvaluator(const UnfrozenEffectiveDepthEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key uf_key_;
  Key depth_key_;
  double alpha_;

 private:
  static Utils::RegisteredFactory<Evaluator, UnfrozenEffectiveDepthEvaluator> fac_;
};

} // namespace Flow
} // namespace Amanzi
