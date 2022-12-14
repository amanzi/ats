/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

//! Evaluates snow melt via USDA - Natural Resources Conservation Service model
/*!

From:  National Engineering Handbook (NEH) part 630, Chapter 11

Uses LandCover for snow_ground_transition parameter.

.. _snow-meltrate-evaluator-spec:
.. admonition:: snow-meltrate-evaluator-spec

   * `"snow melt rate [mm day^-1 C^-1]`" ``[double]`` **2.74**
     the melt rate per degree-day above 0 C.

   * `"air-snow temperature difference [C]`" ``[double]`` **2.0**
     Snow temperature is typicaly a few degrees colder than the air
     temperature at snowmelt. This shifts air temp (positive is colder)
     when calculating the snow temperature.

   * `"surface domain name`" ``[string]`` **SURFACE_DOMAIN** Attempts to
     guess a sane default by the snow domain name.

   KEYS:

   - `"air temperature`"  **SURFACE_DOMAIN-air_temperature**
   - `"snow water equivalent`" **DOMAIN-water_equivalent**


.. note:
    If snow temperature is known, the `"air-snow temperature difference`"
    should be set to 0, and the `"air temperature key`" should be set to
    the snow temperature key instead.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class SnowMeltRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SnowMeltRateEvaluator(Teuchos::ParameterList& plist);
  SnowMeltRateEvaluator(const SnowMeltRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SnowMeltRateEvaluator(*this));
  }

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key temp_key_;
  Key snow_key_;

  double melt_rate_;
  double snow_temp_shift_;

  Key domain_, domain_surf_;
  bool compatibility_checked_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, SnowMeltRateEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
