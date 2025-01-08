/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Richards water content evaluator: the standard form as a function of liquid saturation.
/*!

.. math::
  \Theta = n s \phi V

Specified with evaluator type: `"richards water content`"

.. _field_evaluator_type_richards_water_content-spec:
.. admonition:: field_evaluator_type_richards_water_content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"saturation liquid`"
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_RICHARDS_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_RICHARDS_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class RichardsWaterContentModel;

class RichardsWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit RichardsWaterContentEvaluator(Teuchos::ParameterList& plist);
  RichardsWaterContentEvaluator(const RichardsWaterContentEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;


  Teuchos::RCP<RichardsWaterContentModel> get_model() { return model_; }

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
  Key cv_key_;

  Teuchos::RCP<RichardsWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, RichardsWaterContentEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
