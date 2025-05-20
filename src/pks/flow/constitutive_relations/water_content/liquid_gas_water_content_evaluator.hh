/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Water content for liquid + water vapor.
/*!

.. math::
  \Theta = (n_l s_l + n_g s_g \omega) \phi V


`"evaluator type`" = `"liquid+gas water content`"

.. _evaluator-liquid-gas-water-content-spec:
.. admonition:: evaluator-liquid-gas-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"molar density gas`"
   - `"saturation liquid`"
   - `"saturation gas`"
   - `"mol frac gas`"
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_LIQUID_GAS_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_LIQUID_GAS_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class LiquidGasWaterContentModel;

class LiquidGasWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit LiquidGasWaterContentEvaluator(Teuchos::ParameterList& plist);
  LiquidGasWaterContentEvaluator(const LiquidGasWaterContentEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;


  Teuchos::RCP<LiquidGasWaterContentModel> get_model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

 protected:
  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key sg_key_;
  Key ng_key_;
  Key omega_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidGasWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, LiquidGasWaterContentEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
