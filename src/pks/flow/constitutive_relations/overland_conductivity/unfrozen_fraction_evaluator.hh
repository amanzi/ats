/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the unfrozen fraction model.

*/
/*!

An empirical equation for freezing ponded water -- this is simply a smooth
sinusoidal curve from 0 to 1 over a given transition in temperature.

`"evaluator type`" = `"unfrozen fraction`"

.. _evaluator-unfrozen-fraction-spec:
.. admonition:: evaluator-unfrozen-fraction-spec

   * `"unfrozen fraction model`" ``[unfrozen-fraction-model-spec]``

   KEYS:

   - `"temperature`"
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class UnfrozenFractionModel;

class UnfrozenFractionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  UnfrozenFractionEvaluator(Teuchos::ParameterList& plist);
  UnfrozenFractionEvaluator(const UnfrozenFractionEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  Teuchos::RCP<const UnfrozenFractionModel> get_Model() const { return model_; }
  Teuchos::RCP<UnfrozenFractionModel> get_Model() { return model_; }

 protected:
  Teuchos::RCP<UnfrozenFractionModel> model_;
  Key temp_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, UnfrozenFractionEvaluator> fac_;
};

} // namespace Flow
} // namespace Amanzi

#endif
