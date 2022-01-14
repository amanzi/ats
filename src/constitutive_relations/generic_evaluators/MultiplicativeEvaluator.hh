/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A generic evaluator for multiplying a collection of fields.

/*!

.. _multiplicative-evaluator-spec:
.. admonition:: multiplicative-evaluator-spec
   * `"coefficient`" ``[double]`` **1** A constant prefix to the product.
   * `"enforce positivity`" ``[bool]`` **false** If true, max the result with 0.

   ONE OF
   * `"evaluator dependencies`" ``[Array(string)]`` The fields to multiply.
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
  void Evaluate_(const State& S,
                      const std::vector<CompositeVector*>& result) override;
  void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;

 protected:
  double coef_;
  bool positive_;

 private:
  static Utils::RegisteredFactory<Evaluator,MultiplicativeEvaluator> factory_;
};

} // namespace
} // namespace


