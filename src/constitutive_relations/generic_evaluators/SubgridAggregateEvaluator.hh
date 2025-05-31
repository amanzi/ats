/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

//! SubgridAggregateEvaluator restricts a field to the subgrid version of the same field.
/*!

.. _evaluator-subgrid-aggregate-spec:
.. admonition:: evaluator-subgrid-aggregate-spec

   * `"source domain name`" ``[string]`` Domain name of the source mesh.

   KEYS:

   - `"field`" **SOURCE_DOMAIN-KEY**  Default set from this evaluator's name.

*/


#ifndef AMANZI_RELATIONS_SUBGRID_AGGREGATOR_EVALUATOR_HH_
#define AMANZI_RELATIONS_SUBGRID_AGGREGATOR_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class SubgridAggregateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit SubgridAggregateEvaluator(Teuchos::ParameterList& plist);

  SubgridAggregateEvaluator(const SubgridAggregateEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // custom EnsureEvaluators required to fill dependencies correctly.
  void EnsureEvaluators(State& S) override;

 protected:
  // custom EC required to make sure that this vector is generated on the
  // referencing mesh
  void EnsureCompatibility_Structure_(State& S) override;

  // custom EC required to depend on cells of the subgrid mesh
  void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  void EvaluatePartialDerivative_(const State& S,
                                  const Key& wrt_key,
                                  const Tag& wrt_tag,
                                  const std::vector<CompositeVector*>& result) override;

 protected:
  Key source_domain_;
  Key domain_;
  Key var_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, SubgridAggregateEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
