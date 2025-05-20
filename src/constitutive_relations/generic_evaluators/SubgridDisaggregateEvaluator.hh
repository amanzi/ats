/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! SubgridDisaggregateEvaluator restricts a field to the subgrid version of the same field.
/*!

Note that this evaluator fills exactly one subdomain in a domain set -- there
will be N evaluators each filling one subdomain.

.. _evaluator-subgrid-disaggregate-spec:
.. admonition:: evaluator-subgrid-disaggregate-spec

   * `"source domain name`" ``[string]`` Domain name of the source mesh.

   KEYS:

   - `"field`" **SOURCE_DOMAIN-KEY**  Default set from this evaluator's name.

 */


#ifndef AMANZI_RELATIONS_SUBGRID_DISAGGREGATOR_EVALUATOR_HH_
#define AMANZI_RELATIONS_SUBGRID_DISAGGREGATOR_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class SubgridDisaggregateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit SubgridDisaggregateEvaluator(Teuchos::ParameterList& plist);
  SubgridDisaggregateEvaluator(const SubgridDisaggregateEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
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
  Key domain_index_;
  Key domain_set_;
  Key source_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, SubgridDisaggregateEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
