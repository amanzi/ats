/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Fraction of incoming water that is intercepted.
/*!

Based on CLM 4.5 and Lawrence et al 2007:

.. math::
  I = (P_{rain} + P_{snow}) * \alpha * (1 - exp(-.5(LAI+SAI)))

The interception fraction is everything here after the precip.

.. _interception-fraction-evaluator-spec:
.. admonition:: interception-fraction-evaluator-spec

   * `"interception fraction parameters`" ``[interception-fraction-model-spec]``

   MY KEYS:
   - "interception"
   - "throughfall and drainage rain"
   - "throughfall and drainage snow"

   KEYS:
   - "area index"
   - "precipitation rain"
   - "precipitation snow"
   - "drainage"
   - "air temperature"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class InterceptionFractionModel;

class InterceptionFractionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit InterceptionFractionEvaluator(Teuchos::ParameterList& plist);
  InterceptionFractionEvaluator(const InterceptionFractionEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<InterceptionFractionModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  void InitializeFromPlist_();

 protected:
  Key ai_key_;
  Key rain_key_;
  Key snow_key_;
  Key drainage_key_;
  Key air_temp_key_;

  Key interception_key_;
  Key throughfall_snow_key_;
  Key throughfall_rain_key_;

  Teuchos::RCP<InterceptionFractionModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, InterceptionFractionEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
