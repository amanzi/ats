/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridAggregateEvaluator restricts a field to the subgrid version of the same field.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/


/*!

.. _subgrid-aggregate-evaluator-spec:
.. admonition:: subgrid-aggregate-evaluator-spec

   * `"source domain name`" ``[string]`` Domain name of the source mesh.

   KEYS:
   - `"field`" **SOURCE_DOMAIN-KEY**  Default set from this evaluator's name.

*/


#ifndef AMANZI_RELATIONS_SUBGRID_AGGREGATOR_EVALUATOR_HH_
#define AMANZI_RELATIONS_SUBGRID_AGGREGATOR_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class SubgridAggregateEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  SubgridAggregateEvaluator(Teuchos::ParameterList& plist);

  SubgridAggregateEvaluator(const SubgridAggregateEvaluator& other) = default;
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
  Key var_key_;
  Key source_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SubgridAggregateEvaluator> factory_;
};

} // namespace
} // namespace

#endif

