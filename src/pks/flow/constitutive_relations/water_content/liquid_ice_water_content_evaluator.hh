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
  \Theta = (n_l s_l + n_i s_i) \phi V


Specified with evaluator type: `"liquid+ice water content`"

.. _field_evaluator_type_liquid_ice_water_content-spec:
.. admonition:: field_evaluator_type_liquid_ice_water_content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"molar density ice`"
   - `"saturation liquid`"
   - `"saturation ice`"
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_LIQUID_ICE_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_LIQUID_ICE_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class LiquidIceWaterContentModel;

class LiquidIceWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit LiquidIceWaterContentEvaluator(Teuchos::ParameterList& plist);
  LiquidIceWaterContentEvaluator(const LiquidIceWaterContentEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<LiquidIceWaterContentModel> get_model() { return model_; }

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
  Key si_key_;
  Key ni_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidIceWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, LiquidIceWaterContentEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
