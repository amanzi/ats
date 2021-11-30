/*
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

.. _unfrozen-effective-depth-evaluator-spec
.. admonition:: unfrozen-effective-depth-evaluator-spec

  * `"ice retardation exponent [-]`" ``[double]`` **1.0** exponent alpha
    controlling how quickly ice turns off flow.

  DEPENDENCIES:
  - `"depth`" **DOMAIN-depth**
  - `"unfrozen fraction`" **DOMAIN-unfrozen_fraction**

*/


#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class UnfrozenEffectiveDepthModel;

class UnfrozenEffectiveDepthEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  UnfrozenEffectiveDepthEvaluator(Teuchos::ParameterList& plist);
  UnfrozenEffectiveDepthEvaluator(const UnfrozenEffectiveDepthEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

protected:
  Key uf_key_;
  Key depth_key_;
  double alpha_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,UnfrozenEffectiveDepthEvaluator> fac_;


};

} //namespace
} //namespace


