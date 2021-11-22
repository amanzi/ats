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
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class MultiplicativeEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  MultiplicativeEvaluator(Teuchos::ParameterList& plist);
  MultiplicativeEvaluator(const MultiplicativeEvaluator& other) = default;

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  void EvaluateField_(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<CompositeVector>& result);
  void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  double coef_;
  bool positive_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,MultiplicativeEvaluator> factory_;
};

} // namespace
} // namespace


