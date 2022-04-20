/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridDisaggregateEvaluator restricts a field to the subgrid version of the same field.

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

.. _subgrid-disaggregate-evaluator-spec:
.. admonition:: subgrid-disaggregate-evaluator-spec

   * `"source domain name`" ``[string]`` Domain name of the source mesh.

   KEYS:
   - `"field`" **SOURCE_DOMAIN-KEY**  Default set from this evaluator's name.


 */


#ifndef AMANZI_RELATIONS_SUBGRID_DISAGGREGATOR_EVALUATOR_HH_
#define AMANZI_RELATIONS_SUBGRID_DISAGGREGATOR_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class SubgridDisaggregateEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  SubgridDisaggregateEvaluator(Teuchos::ParameterList& plist);

  SubgridDisaggregateEvaluator(const SubgridDisaggregateEvaluator& other) = default;
  Teuchos::RCP<FieldEvaluator> Clone() const;

  void
  EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:

  // Required methods from SecondaryVariableFieldEvaluator
  void EvaluateField_(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<CompositeVector>& result);
  void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);



 protected:
  Key source_domain_;
  Key domain_;
  Key source_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SubgridDisaggregateEvaluator> factory_;
};

} // namespace
} // namespace

#endif

