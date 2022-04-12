/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! SubgridMobileDepthEvaluator: calculates mobile depth including a depression storage term.
/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Effectively, this is

.. math:
   \delta_{mobile} = max(0, \delta - \delta_{depression})

from Jan et al WRR 2018.

It is a bit trickier than this because depression depth is only provided on
cells, while mobile depth (and ponded depth) are on both cells and boundary
faces.

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class SubgridMobileDepthEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SubgridMobileDepthEvaluator(Teuchos::ParameterList& plist);
  SubgridMobileDepthEvaluator(const SubgridMobileDepthEvaluator& other) = default;

  Teuchos::RCP<FieldEvaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

private:
  Key depth_key_;
  Key depr_depth_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SubgridMobileDepthEvaluator> factory_;
};

} //namespace
} //namespace
} //namespace


