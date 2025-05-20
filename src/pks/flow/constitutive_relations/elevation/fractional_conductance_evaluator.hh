/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan
*/

//! FractionalConductanceEvaluator: an obstruction-drag factor.
/*!

This implements the term,

.. math:
    \frac{\Phi(\delta) - \Phi(\delta_d)}{\delta - \delta_d}

from Jan et al WRR 2018.

.. _evaluator-fractional-conductance-spec
.. admonition:: evaluator-fractional-conductance-spec

  DEPENDENCIES:
  - `"volumetric ponded depth`" **DOMAIN-volumetric_ponded_depth**
  - `"mobile depth`" **DOMAIN-mobile_depth**  Note, this is typically d - d_{depr}
  - `"microtopographic relief`" **DOMAIN-microtopographic_relief**
  - `"excluded volume`" **DOMAIN-excluded_volume**
  - `"depression depth`" **DOMAIN-depression_depth**

*/

#ifndef AMANZI_FLOWRELATIONS_FRACTIONAL_CONDUCTANCE_EVALUATOR_
#define AMANZI_FLOWRELATIONS_FRACTIONAL_CONDUCTANCE_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class FractionalConductanceEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit FractionalConductanceEvaluator(Teuchos::ParameterList& plist);
  FractionalConductanceEvaluator(const FractionalConductanceEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 private:
  Key mobile_depth_key_;
  Key vpd_key_;
  Key depr_depth_key_;
  Key delta_ex_key_, delta_max_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, FractionalConductanceEvaluator> factory_;
};

} // namespace FlowRelations
} // namespace Flow
} // namespace Amanzi

#endif
