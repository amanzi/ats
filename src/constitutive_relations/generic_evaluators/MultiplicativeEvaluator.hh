/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A generic evaluator for multiplying a collection of fields.
/*!

.. _evaluator-multiplicative-spec:
.. admonition:: evaluator-multiplicative-spec
   * `"coefficient`" ``[double]`` **1** A constant prefix to the product.

   * `"enforce positivity`" ``[bool]`` **false** If true, max the result with 0.

   * `"DEPENDENCY dof`" ``[double]`` **0** Degree of Freedom for each given
     dependency to use in the multiplication.  NOTE, this should only be
     provided if the dependency has more than 1 DoF -- if it just has one
     leaving this blank results in better error checking than providing the
     value 0 manually.

   ONE OF

   * `"dependencies`" ``[Array(string)]`` The fields to multiply.

   OR

   * `"evaluator dependency suffixes`" ``[Array(string)]``

   END

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class MultiplicativeEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit MultiplicativeEvaluator(Teuchos::ParameterList& plist);
  MultiplicativeEvaluator(const MultiplicativeEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  void EvaluatePartialDerivative_(const State& S,
                                  const Key& wrt_key,
                                  const Tag& wrt_tag,
                                  const std::vector<CompositeVector*>& result) override;

  void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  double coef_;
  bool positive_;

  std::vector<int> dofs_;
  std::vector<bool> dof_provided_;
  bool any_dof_provided_;

 private:
  static Utils::RegisteredFactory<Evaluator, MultiplicativeEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi
