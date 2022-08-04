/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
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

.. _fractional-conductance-evaluator-spec
.. admonition:: fractional-conductance-evaluator-spec

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
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class FractionalConductanceEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  FractionalConductanceEvaluator(Teuchos::ParameterList& plist);
  FractionalConductanceEvaluator(const FractionalConductanceEvaluator& other) = default;

  Teuchos::RCP<FieldEvaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

private:
  Key mobile_depth_key_;
  Key vpd_key_;
  Key depr_depth_key_;
  Key delta_ex_key_, delta_max_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,FractionalConductanceEvaluator> factory_;
};

} //namespace
} //namespace
} //namespace

#endif
