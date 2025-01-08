/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Three phase water content: vapor, liquid, and ice.
/*!

.. math::
  \Theta = (n_l s_l + n_i s_i + n_g s_g \omega_g ) \phi V

Specified with evaluator type: `"three phase water content`"

.. _field_evaluator_type_three_phase_water_content-spec:
.. admonition:: field_evaluator_type_three_phase_water_content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"saturation liquid`"
   - `"molar density ice`"
   - `"saturation ice`"
   - `"molar density gas`"
   - `"saturation gas`"
   - `"molar fraction gas`"
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ThreePhaseWaterContentModel;

class ThreePhaseWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ThreePhaseWaterContentEvaluator(Teuchos::ParameterList& plist);
  ThreePhaseWaterContentEvaluator(const ThreePhaseWaterContentEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;


  Teuchos::RCP<ThreePhaseWaterContentModel> get_model() { return model_; }

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
  Key sg_key_;
  Key ng_key_;
  Key omega_key_;
  Key cv_key_;

  Teuchos::RCP<ThreePhaseWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, ThreePhaseWaterContentEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
