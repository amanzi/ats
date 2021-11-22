/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A generic evaluator for summing a collection of fields.

/*!

.. _additive-evaluator-spec:
.. admonition:: additive-evaluator-spec
   * `"constant shift`" ``[double]`` **0** A constant value to add to the sum.

   * `"DEPENDENCY coefficient`" ``[double]`` A multiple for each dependency in
     the list below.

   ONE OF
   * `"evaluator dependencies`" ``[Array(string)]`` The things to sum.
   OR
   * `"evaluator dependency suffixes`" ``[Array(string)]``
   END

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class AdditiveEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  AdditiveEvaluator(Teuchos::ParameterList& plist);

  AdditiveEvaluator(const AdditiveEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  void EvaluateField_(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<CompositeVector>& result);
  void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  std::map<Key, double> coefs_;
  double shift_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,AdditiveEvaluator> factory_;
};

} // namespace
} // namespace


