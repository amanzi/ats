/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Partitions water between interception and throughfall, combining throughfall with drainage.
/*!

Based on CLM 4.5 and Lawrence et al 2007, interception is given by:

.. math::
   I = (P_{rain} + P_{snow}) * \alpha * (1 - exp(-.5(LAI)))

Throughfall is given by:

.. math::
   T = (P_{rain} + P_{snow}) - I

Drainage is provided as input here, as a total drainage from the canopy.  The
phase of this drainage is assumed to match the phase of the precipitation.  So
if it is raining, drainage is rain, while if it is 50/50 rain and snow,
drainage is also 50/50 rain and snow.  If total precipitation is 0, then
drainage is partitioned by air temperature (above 0C --> all rain, otherwise
all snow).  This evaluator partitions the drainage and sums it with throughfall
to compute the total source, in each phase, to the layer below the canopy (snow
and/or ground surface).

.. _interception-fraction-evaluator-spec:
.. admonition:: interception-fraction-evaluator-spec

   * `"interception fraction parameters`" ``[interception-fraction-model-spec]``

   MY KEYS:
   - `"interception`" **DOMAIN-interception**
   - `"throughfall and drainage rain`" **DOMAIN-throughfall_drainage_rain**
   - `"throughfall and drainage snow`" **DOMAIN-throughfall_drainage_snow**

   KEYS:
   - `"area index`" **DOMAIN-area_index**
   - `"precipitation rain`" **DOMAIN_SURFACE-precipitation_rain**
   - `"precipitation snow`" **DOMAIN_SNOW-precipitation**
   - `"drainage`" **DOMAIN-drainage**
   - `"air temperature`" **DOMAIN_SURFACE-air_temperature**

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
