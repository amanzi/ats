/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ahmad Jan (jana@ornl.gov)
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
#include "secondary_variable_field_evaluator.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class SnowMeltRateEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SnowMeltRateEvaluator(Teuchos::ParameterList& plist);
  SnowMeltRateEvaluator(const SnowMeltRateEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new SnowMeltRateEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

 protected:

  Key temp_key_;
  Key snow_key_;

  double melt_rate_;
  double snow_temp_shift_;

  Key domain_, domain_surf_;
  bool compatibility_checked_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SnowMeltRateEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

