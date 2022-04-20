/*
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
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class InterceptionFractionModel;

class InterceptionFractionEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  explicit
  InterceptionFractionEvaluator(Teuchos::ParameterList& plist);
  InterceptionFractionEvaluator(const InterceptionFractionEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) override;

  Teuchos::RCP<InterceptionFractionModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

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
  static Utils::RegisteredFactory<FieldEvaluator,InterceptionFractionEvaluator> reg_;
};

} //namespace
} //namespace
} //namespace


